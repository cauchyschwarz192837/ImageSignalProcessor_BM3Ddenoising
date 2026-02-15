#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

static float sigma = 25.0f;
static float lamb2d = 2.0f;
static float lamb3d = 2.7f;

static int Step1_ThreDist = 2500;
static int Step1_MaxMatch = 16;
static int Step1_BlockSize = 8;
static int Step1_spdup_factor = 3;
static int Step1_WindowSize = 39;

static int Step2_ThreDist = 400;
static int Step2_MaxMatch = 32;
static int Step2_BlockSize = 8;
static int Step2_spdup_factor = 3;
static int Step2_WindowSize = 39;

static float Kaiser_Window_beta = 2.0f;

// Helper funcs

static double I0(double x) {
    double ax = std::abs(x);
    if (ax < 3.75) {
        double t = x / 3.75;
        double t2 = t * t;
        return 1.0 + t2 * (3.5156229 + t2*(3.0899424 + t2*(1.2067492 + t2 * (0.2659732 + t2*(0.0360768 + t2*0.0045813)))));
    } else {
        double t = 3.75 / ax;
        return (std::exp(ax) / std::sqrt(ax)) * (0.39894228 + t*(0.01328592 + t * (0.00225319 + t*(-0.00157565 + t*(0.00916281 + t*(-0.02057706 + t * (0.02635537 + t*(-0.01647633 + t*0.00392377))))))));
    }
}

static std::vector<float> kaiser1D(int N, float beta) {
    std::vector<float> w(N);
    if (N == 1) { 
        w[0] = 1.0f; 
        return w; 
    }
    double denom = I0(beta);
    for (int n = 0; n < N; n++) {
        double r = 2.0 * n / (N - 1) - 1.0;  // in [-1, 1]
        double val = I0(beta * std::sqrt(1.0 - r * r)) / denom;
        w[n] = static_cast<float>(val);
    }
    return w;
}

static cv::Mat kaiser2D(int B, float beta) {
    auto w = kaiser1D(B, beta);
    cv::Mat K(B, B, CV_32F);
    for (int i = 0; i < B; i++) {
        for (int j = 0; j < B; j++) {
            K.at<float>(i, j) = w[i] * w[j];
        }
    }
    return K;
}

static cv::Mat addNoiseGaussian(const cv::Mat& img_u8, float sigma_std) {
    cv::Mat imgf; img_u8.convertTo(imgf, CV_32F);
    cv::Mat noise(imgf.size(), CV_32F);
    cv::randn(noise, 0.0, sigma_std);
    return imgf + noise;
}

static double computePSNR_u8_vs_float(const cv::Mat& img_u8, const cv::Mat& img_float) {
    CV_Assert(img_u8.size() == img_float.size());
    cv::Mat a; img_u8.convertTo(a, CV_32F);
    cv::Mat diff = a - img_float;
    cv::Mat diff2; cv::multiply(diff, diff, diff2);
    double mse = cv::sum(diff2)[0] / (double)img_u8.total();
    if (mse <= 1e-12) {
        return 99.0;
    }
    double rmse = std::sqrt(mse);
    return 20.0 * std::log10(255.0 / rmse);
}

static inline cv::Rect clampRect(int x, int y, int w, int h, int W, int H) {
    x = std::max(0, std::min(x, W - w));
    y = std::max(0, std::min(y, H - h));
    return cv::Rect(x, y, w, h);
}

static cv::Point searchWindowTopLeft(const cv::Mat& img, const cv::Point& ref, int blockSize, int windowSize) {
    if (blockSize >= windowSize) {
        throw std::runtime_error("BlockSize must be smaller than WindowSize");
    }

    int H = img.rows, W = img.cols;
    int sx = std::max(0, ref.x + (blockSize - windowSize) / 2);
    int sy = std::max(0, ref.y + (blockSize - windowSize) / 2);

    // ensure window fits
    sx = std::min(sx, H - windowSize);
    sy = std::min(sy, W - windowSize);

    sx = std::max(0, sx);
    sy = std::max(0, sy);
    return cv::Point(sx, sy);
}

static cv::Mat dct2D(const cv::Mat& block32f) {
    cv::Mat out;
    cv::dct(block32f, out, cv::DCT_ROWS);
    cv::dct(block32f, out);
    return out;
}

static cv::Mat idct2D(const cv::Mat& block32f) {
    cv::Mat out;
    cv::idct(block32f, out);
    return out;
}

struct PreDCTTable {
    int H = 0, W = 0, B = 0;
    int nH = 0, nW = 0;
    std::vector<cv::Mat> blocks; // each is BxB CV_32F

    inline int idx(int i, int j) const { 
        return i * nW + j; 
    }

    const cv::Mat& at(int i, int j) const { 
        return blocks[idx(i, j)]; 
    }

    cv::Mat& at(int i, int j) { return blocks[idx(i, j)]; }
};

static PreDCTTable preDCT(const cv::Mat& img32f, int B) {
    PreDCTTable T;
    T.H = img32f.rows; T.W = img32f.cols; T.B = B;
    T.nH = T.H - B + 1;
    T.nW = T.W - B + 1;
    if (T.nH <= 0 || T.nW <= 0) {
        throw std::runtime_error("Image too small for given BlockSize");
    }

    T.blocks.resize((size_t)T.nH * (size_t)T.nW);
    for (int i = 0; i < T.nH; i++) {
        for (int j = 0; j < T.nW; j++) {
            cv::Mat block = img32f(cv::Rect(j, i, B, B)); // note: Rect(x,y,w,h) => (col,row)
            T.at(i, j) = dct2D(block);
        }
    }
    return T;
}

// (1)

static float step1_computeDist(cv::Mat A, cv::Mat B) {
    // thresholding for sigma>40
    if (sigma > 40.0f) {
        float thre = lamb2d * sigma;
        cv::Mat maskA = cv::abs(A) < thre;
        cv::Mat maskB = cv::abs(B) < thre;
        A.setTo(0.0f, maskA);
        B.setTo(0.0f, maskB);
    }

    cv::Mat diff = A - B;
    double nrm = cv::norm(diff, cv::NORM_L2);
    int BS = A.rows;
    return (float)((nrm * nrm) / (double)(BS * BS));
}

struct Group1 {
    std::vector<cv::Point> pos; // (i,j) in DCT-table coords (top-left block location)
    std::vector<cv::Mat>   grp; // DCT blocks
};

static Group1 step1_grouping(const cv::Mat& noisy32f, const cv::Point& ref, const PreDCTTable& dctAll, int B, int threDist, int maxMatch, int windowSize) {
    cv::Point winTL = searchWindowTopLeft(noisy32f, ref, B, windowSize);
    int searchSpan = windowSize - B + 1;
    const cv::Mat& refD = dctAll.at(ref.x, ref.y);

    struct Cand { float dist; cv::Point p; };
    std::vector<Cand> cands;
    cands.reserve((size_t)searchSpan * (size_t)searchSpan);

    for (int di = 0; di < searchSpan; di++) {
        for (int dj = 0; dj < searchSpan; dj++) {
            int bi = winTL.x + di;
            int bj = winTL.y + dj;
            const cv::Mat& sD = dctAll.at(bi, bj);
            float dist = step1_computeDist(refD, sD);
            if (dist < (float)threDist) {
                cands.push_back({dist, cv::Point(bi, bj)});
            }
        }
    }

    if ((int)cands.size() > maxMatch) {
        std::nth_element(cands.begin(), cands.begin() + maxMatch, cands.end(), [](const Cand& a, const Cand& b){ return a.dist < b.dist; });
        cands.resize(maxMatch);
    }

    Group1 out;
    out.pos.reserve(cands.size());
    out.grp.reserve(cands.size());
    for (auto& c : cands) {
        out.pos.push_back(c.p);
        out.grp.push_back(dctAll.at(c.p.x, c.p.y).clone());
    }
    return out;
}

static int step1_3DFiltering(std::vector<cv::Mat>& groupDCT) {
    int N = (int)groupDCT.size();
    if (N == 0) {
        return 0;
    }

    int B = groupDCT[0].rows;

    int L = (N % 2 == 0) ? N : (N + 1);

    float thre = lamb3d * sigma;
    int nonzero_cnt = 0;

    cv::Mat vec(L, 1, CV_32F, cv::Scalar(0));
    cv::Mat vecD, vecID;

    for (int i = 0; i < B; i++) {
        for (int j = 0; j < B; j++) {
            // fill first N, pad last with 0 if needed
            for (int k = 0; k < N; ++k) {
                vec.at<float>(k, 0) = groupDCT[k].at<float>(i, j);
            }
            if (L > N) {
                vec.at<float>(N, 0) = 0.0f;
            }

            cv::dct(vec, vecD);

            // hard threshold only on first N meaningful coeffs
            for (int k = 0; k < N; k++) {
                float& v = vecD.at<float>(k, 0);
                if (std::abs(v) < thre) {
                    v = 0.0f;
                } else {
                    nonzero_cnt++;
                }
            }
            if (L > N) {
                vecD.at<float>(N, 0) = 0.0f;
            }

            cv::idct(vecD, vecID);

            // write back only first N
            for (int k = 0; k < N; ++k) {
                groupDCT[k].at<float>(i, j) = vecID.at<float>(k, 0);
            }
        }
    }
    return nonzero_cnt;
}

static void step1_aggregate(const std::vector<cv::Mat>& groupDCT, const std::vector<cv::Point>& positions, cv::Mat& basicImg, cv::Mat& basicW, const cv::Mat& kaiser, int nonzero_cnt) {
    if (groupDCT.empty()) {
        return;
    }

    float baseWeight;
    if (nonzero_cnt < 1) {
        baseWeight = 1.0f;
    } else {
        baseWeight = 1.0f / (sigma * sigma * (float)nonzero_cnt);
    }

    cv::Mat blockWeight = baseWeight * kaiser; // B x B

    int B = groupDCT[0].rows;

    for (size_t t = 0; t < positions.size(); t++) {
        cv::Mat spatial = idct2D(groupDCT[t]);

        int bi = positions[t].x;
        int bj = positions[t].y;

        cv::Rect roi(bj, bi, B, B);
        basicImg(roi) += blockWeight.mul(spatial);
        basicW(roi) += blockWeight;
    }
}

static cv::Mat bm3d_step1(const cv::Mat& noisy32f) {
    int B = Step1_BlockSize;

    cv::Mat basicImg = cv::Mat::zeros(noisy32f.size(), CV_32F);
    cv::Mat basicW = cv::Mat::zeros(noisy32f.size(), CV_32F);
    cv::Mat kaiser = kaiser2D(B, Kaiser_Window_beta);

    PreDCTTable dctAll = preDCT(noisy32f, B);

    for (int ii = 0; ii < (int)((dctAll.nH - 1) / Step1_spdup_factor) + 2; ii++) {
        for (int jj = 0; jj < (int)((dctAll.nW - 1) / Step1_spdup_factor) + 2; jj++) {
            int ri = std::min(Step1_spdup_factor * ii, dctAll.nH - 1);
            int rj = std::min(Step1_spdup_factor * jj, dctAll.nW - 1);
            cv::Point ref(ri, rj);

            Group1 g = step1_grouping(noisy32f, ref, dctAll, B, Step1_ThreDist,
                                      Step1_MaxMatch, Step1_WindowSize);

            int nonzero = step1_3DFiltering(g.grp);
            step1_aggregate(g.grp, g.pos, basicImg, basicW, kaiser, nonzero);
        }
    }

    // avoid div by zero
    cv::Mat wSafe = basicW.clone();
    wSafe.setTo(1.0f, wSafe == 0.0f);

    cv::Mat out = basicImg / wSafe;
    return out;
}

// (2)

static float step2_computeDist(const cv::Mat& basic32f, const cv::Point& p1, const cv::Point& p2, int B) {
    cv::Rect r1(p1.y, p1.x, B, B);
    cv::Rect r2(p2.y, p2.x, B, B);
    cv::Mat a = basic32f(r1);
    cv::Mat b = basic32f(r2);
    cv::Mat diff = a - b;
    double nrm = cv::norm(diff, cv::NORM_L2);
    return (float)((nrm * nrm) / (double)(B * B));
}

struct Group2 {
    std::vector<cv::Point> pos;
    std::vector<cv::Mat> grpBasicDCT;
    std::vector<cv::Mat> grpNoisyDCT;
};

static Group2 step2_grouping(const cv::Mat& basic32f, const cv::Mat& noisy32f, const cv::Point& ref, int B, int threDist, int maxMatch, int windowSize, const PreDCTTable& dctBasic, const PreDCTTable& dctNoisy) {
    cv::Point winTL = searchWindowTopLeft(basic32f, ref, B, windowSize);
    int searchSpan = windowSize - B + 1;

    struct Cand { 
        float dist; cv::Point p; 
    };

    std::vector<Cand> cands;
    cands.reserve((size_t)searchSpan * (size_t)searchSpan);

    for (int di = 0; di < searchSpan; di++) {
        for (int dj = 0; dj < searchSpan; dj++) {
            cv::Point p(winTL.x + di, winTL.y + dj);
            float dist = step2_computeDist(basic32f, ref, p, B);
            if (dist < (float)threDist) {
                cands.push_back({dist, p});
            }
        }
    }

    if ((int)cands.size() > maxMatch) {
        std::nth_element(cands.begin(), cands.begin() + maxMatch, cands.end(), [](const Cand& a, const Cand& b){ return a.dist < b.dist; });
        cands.resize(maxMatch);
    }

    Group2 out;
    out.pos.reserve(cands.size());
    out.grpBasicDCT.reserve(cands.size());
    out.grpNoisyDCT.reserve(cands.size());

    for (auto& c : cands) {
        out.pos.push_back(c.p);
        out.grpBasicDCT.push_back(dctBasic.at(c.p.x, c.p.y).clone());
        out.grpNoisyDCT.push_back(dctNoisy.at(c.p.x, c.p.y).clone());
    }
    return out;
}

static float step2_3DFiltering(std::vector<cv::Mat>& grpBasicDCT, std::vector<cv::Mat>& grpNoisyDCT) {
    int N = (int)grpNoisyDCT.size();
    if (N == 0) {
        return 1.0f;
    }
    int B = grpNoisyDCT[0].rows;

    int L = (N % 2 == 0) ? N : (N + 1);

    float WeightSum = 0.0f;
    float coef = 1.0f / (float)N;

    cv::Mat vb(L, 1, CV_32F, cv::Scalar(0));
    cv::Mat vn(L, 1, CV_32F, cv::Scalar(0));
    cv::Mat vbD, vnD, vnID;

    for (int i = 0; i < B; i++) {
        for (int j = 0; j < B; j++) {
            for (int k = 0; k < N; k++) {
                vb.at<float>(k, 0) = grpBasicDCT[k].at<float>(i, j);
                vn.at<float>(k, 0) = grpNoisyDCT[k].at<float>(i, j);
            }

            if (L > N) { 
                vb.at<float>(N, 0) = 0.0f; vn.at<float>(N, 0) = 0.0f; 
            }

            cv::dct(vb, vbD);
            cv::dct(vn, vnD);

            for (int k = 0; k < N; k++) {
                float b = vbD.at<float>(k, 0);
                float val = (b*b) * coef;
                val = val / (val + sigma*sigma);
                vnD.at<float>(k, 0) *= val;
                WeightSum += val;
            }

            if (L > N) {
                vnD.at<float>(N, 0) = 0.0f;
            }

            cv::idct(vnD, vnID);

            for (int k = 0; k < N; ++k) {
                grpNoisyDCT[k].at<float>(i, j) = vnID.at<float>(k, 0);
            }
        }
    }

    return (WeightSum > 0.0f) ? (1.0f / (sigma*sigma * WeightSum)) : 1.0f;
}

static void step2_aggregate(const std::vector<cv::Mat>& grpNoisyDCT, const std::vector<cv::Point>& positions, float wienerWeight, cv::Mat& finalImg, cv::Mat& finalW, const cv::Mat& kaiser) {
    if (grpNoisyDCT.empty()) {
        return;
    }
    int B = grpNoisyDCT[0].rows;

    cv::Mat blockWeight = wienerWeight * kaiser;

    for (size_t t = 0; t < positions.size(); t++) {
        cv::Mat spatial = idct2D(grpNoisyDCT[t]);

        int bi = positions[t].x;
        int bj = positions[t].y;

        cv::Rect roi(bj, bi, B, B);
        finalImg(roi) += blockWeight.mul(spatial);
        finalW(roi) += blockWeight;
    }
}

static cv::Mat bm3d_step2(const cv::Mat& basic32f, const cv::Mat& noisy32f) {
    int B = Step2_BlockSize;

    cv::Mat finalImg = cv::Mat::zeros(basic32f.size(), CV_32F);
    cv::Mat finalW = cv::Mat::zeros(basic32f.size(), CV_32F);
    cv::Mat kaiser = kaiser2D(B, Kaiser_Window_beta);

    PreDCTTable dctNoisy = preDCT(noisy32f, B);
    PreDCTTable dctBasic = preDCT(basic32f, B);

    for (int ii = 0; ii < (int)((dctBasic.nH - 1) / Step2_spdup_factor) + 2; ii++) {
        for (int jj = 0; jj < (int)((dctBasic.nW - 1) / Step2_spdup_factor) + 2; jj++) {
            int ri = std::min(Step2_spdup_factor * ii, dctBasic.nH - 1);
            int rj = std::min(Step2_spdup_factor * jj, dctBasic.nW - 1);
            cv::Point ref(ri, rj);

            Group2 g = step2_grouping(basic32f, noisy32f, ref, B, Step2_ThreDist, Step2_MaxMatch, Step2_WindowSize, dctBasic, dctNoisy);

            float wienerWeight = step2_3DFiltering(g.grpBasicDCT, g.grpNoisyDCT);
            step2_aggregate(g.grpNoisyDCT, g.pos, wienerWeight, finalImg, finalW, kaiser);
        }
    }

    cv::Mat wSafe = finalW.clone();
    wSafe.setTo(1.0f, wSafe == 0.0f);

    return finalImg / wSafe;
}

//-------------------------------------------------

int main(int argc, char** argv) {
    std::string inPath = "billevans.png";
    if (argc >= 2) {
        inPath = argv[1];
    }

    cv::Mat img = cv::imread(inPath, cv::IMREAD_GRAYSCALE);
    if (img.empty()) {
        std::cerr << "Failed to read image: " << inPath << "\n";
        return 1;
    }

    cv::Mat noisy = addNoiseGaussian(img, sigma);
    double t0 = (double)cv::getTickCount();

    cv::Mat basic = bm3d_step1(noisy);
    double basicPSNR = computePSNR_u8_vs_float(img, basic);
    std::cout << "Basic PSNR: " << basicPSNR << " dB\n";

    cv::Mat basic_u8;
    cv::normalize(basic, basic_u8, 0, 255, cv::NORM_MINMAX);
    basic_u8.convertTo(basic_u8, CV_8U);
    cv::imwrite("basicBill.png", basic_u8);

    double t1 = (double)cv::getTickCount();
    std::cout << "Step1 time (s): " << (t1 - t0) / cv::getTickFrequency() << "\n";

    cv::Mat final = bm3d_step2(basic, noisy);
    double finalPSNR = computePSNR_u8_vs_float(img, final);
    std::cout << "Final PSNR: " << finalPSNR << " dB\n";

    cv::Mat final_u8;
    cv::normalize(final, final_u8, 0, 255, cv::NORM_MINMAX);
    final_u8.convertTo(final_u8, CV_8U);
    cv::imwrite("finalBill.png", final_u8);

    double t2 = (double)cv::getTickCount();
    std::cout << "Step2 time (s): " << (t2 - t1) / cv::getTickFrequency() << "\n";

    return 0;
}