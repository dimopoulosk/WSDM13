#ifndef _PFOR_H_
#define _PFOR_H_

typedef void (*pf)(unsigned int *p, unsigned int *w);

void pack(unsigned int *v, unsigned int b, unsigned int n, unsigned int *w);
int pack_encode(unsigned int **w, unsigned int *p, int num);

unsigned int *pack_decode(unsigned int *_p, unsigned int *_w, int flag);
unsigned int *pack_decode(unsigned int *_p, unsigned int *_w, int flag,int base);

void pack(unsigned int *v, unsigned int b, unsigned int n, unsigned int *w);

void unpack2(unsigned int *p, unsigned int *w);
void unpack3(unsigned int *p, unsigned int *w);
void unpack4(unsigned int *p, unsigned int *w);
void unpack5(unsigned int *p, unsigned int *w);
void unpack6(unsigned int *p, unsigned int *w);
void unpack7(unsigned int *p, unsigned int *w);
void unpack8(unsigned int *p, unsigned int *w);
void unpack9(unsigned int *p, unsigned int *w);
void unpack10(unsigned int *p, unsigned int *w);
void unpack12(unsigned int *p, unsigned int *w);
void unpack16(unsigned int *p, unsigned int *w);
void unpack20(unsigned int *p, unsigned int *w);
void unpack32(unsigned int *p, unsigned int *w);

void prefix(unsigned int *_p, unsigned int s);
void prefix2(unsigned int *p);
void prefix3(unsigned int *p);

#endif
