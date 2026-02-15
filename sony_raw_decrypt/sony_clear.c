#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <netinet/in.h>

#ifdef htonl
#undef htonl
#endif

#ifdef ntohl
#undef ntohl
#endif

typedef unsigned char uchar;

typedef union {
	unsigned char uc[4];
	unsigned long ul;
} UDWord;

unsigned long htonl(unsigned long ll) {
	UDWord x, y;

	x.ul = ll;
	y.uc[0] = x.uc[3];
	y.uc[1] = x.uc[2];
	y.uc[2] = x.uc[1];
	y.uc[3] = x.uc[0];
	return y.ul;
}

unsigned long ntohl(unsigned long ll) {
	return htonl(ll);
}

void sony_decrypt (void *buf, int len, int key) {
  unsigned pad[128], *data = buf;
  int i;

  for (i = 0; i < 4; i++) {
    pad[i] = key = key * 48828125 + 1;
  }
  pad[3] = pad[3] << 1 | (pad[0]^pad[2]) >> 31;
  for (i = 4; i < 127; i++) {
    pad[i] = (pad[i-4]^pad[i-2]) << 1 | (pad[i-3]^pad[i-1]) >> 31;
  }
  for (i = 0; i < 127; i++) {
    pad[i] = htonl(pad[i]);
  }
  for (; i < len+127; i++, data++) {
    *data ^= pad[i & 127] = pad[(i + 1) & 127] ^ pad[(i + 65) & 127];
  }
}

void sony_clear (uchar *buffer, int length) {
  unsigned *ip, key0, key1 = 0, key2 = 0;
  uchar *cp;

  cp = buffer + 200896;
  ip = (void *) cp;
  key0 = ntohl(ip[*cp]);
  sony_decrypt (buffer + 164600, 9074, key0);
  for (int i = 4; i--;) {
    key1 = key1 << 8 | buffer[164610 + i];
    key2 = key2 << 8 | buffer[164622 + i];
  }
  sony_decrypt (buffer + 164640, 174376, key1);
  sony_decrypt (buffer + 862144, (length - 862144) / 4, key2);
}

int main (int argc, char **argv) {
  FILE *fp;
  char name[512], *buffer;
  int length;

  for (int arg = 1; arg < argc; arg++) {
    fp = fopen(argv[arg], "rb");

    if (!fp) {
      perror(argv[arg]);
      continue;
    }

    fseek(fp, 0, SEEK_END);  // go to end of file
    length = ftell(fp);  // total file size in bytes
    if (length < 0x100000) {  // smaller than 1 MB
      fprintf(stderr, "%s is too small!\n", argv[arg]);
      fclose(fp);
      continue;
    }

    buffer = malloc(length);
    if (!buffer) {
      fprintf (stderr, "%s is too big!\n", argv[arg]);
      fclose(fp);
      continue;
    }

    // read into RAM buffer
    fseek(fp, 0, SEEK_SET);
    // size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
    fread(buffer, 1, length, fp);  // 1 byte at a time
    fclose(fp);

    // decrypt
    sony_clear(buffer, length);

    strcpy(name, argv[arg]);
    strcat(name, ".clear");

    fp = fopen(name, "wb");
    fwrite(buffer, 1, length, fp);  // 1 byte at a time
    free(buffer); // free memory
    fclose(fp);
  }

  return 0;
}
