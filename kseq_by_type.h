/* The MIT License

   Copyright (c) 2008, 2009, 2011 Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Last Modified: 05MAR2012 */

//#ifndef AC_KSEQ_H
//#define AC_KSEQ_H

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#define KS_SEP_MAX   2

#define GLUE_HELPER(x, y) x##y
#define GLUE(x, y) GLUE_HELPER(x, y)

#define __KS_TYPE(type_t, reader)						\
	typedef struct GLUE(__kstream_t_, reader) {				\
		unsigned char *buf;						\
		int begin, end, is_eof;					\
		type_t f;								\
	} GLUE(kstream_t_, reader);

#define ks_err(ks) ((ks)->end == -1)
#define ks_eof(ks) ((ks)->is_eof && (ks)->begin >= (ks)->end)
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)

#define __KS_BASIC(type_t, __bufsize, reader)								\
	GLUE(kstream_t_, reader) *GLUE(ks_init_, reader)(type_t f)						\
	{																\
		GLUE(kstream_t_, reader) *ks = (GLUE(kstream_t_, reader)*)calloc(1, sizeof(GLUE(kstream_t_, reader)));	\
		ks->f = f;													\
		ks->buf = (unsigned char*)malloc(__bufsize);				\
		return ks;													\
	}																\
	void GLUE(ks_destroy_, reader)(GLUE(kstream_t_, reader) *ks)					\
	{																\
		if (ks) {													\
			free(ks->buf);											\
			free(ks);												\
		}															\
	}

#define __KS_GETC(__read, __bufsize, reader)						\
	int GLUE(ks_getc_, reader)(GLUE(kstream_t_, reader) *ks)				\
	{														\
		if (ks_err(ks)) return -3;							\
		if (ks->is_eof && ks->begin >= ks->end) return -1;	\
		if (ks->begin >= ks->end) {							\
			ks->begin = 0;									\
			ks->end = __read(ks->f, ks->buf, __bufsize);	\
			if (ks->end == 0) { ks->is_eof = 1; return -1;}	\
			if (ks->end == -1) { ks->is_eof = 1; return -3;}\
		}													\
		return (int)ks->buf[ks->begin++];					\
	}

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#define __KS_GETUNTIL(__read, __bufsize, reader)								\
	int GLUE(ks_getuntil2_, reader)(GLUE(kstream_t_, reader) *ks, int delimiter, kstring_t *str, int *dret, int append) \
	{																	\
		int gotany = 0;													\
		if (dret) *dret = 0;											\
		str->l = append? str->l : 0;									\
		for (;;) {														\
			int i;														\
			if (ks_err(ks)) return -3;									\
			if (ks->begin >= ks->end) {									\
				if (!ks->is_eof) {										\
					ks->begin = 0;										\
					ks->end = __read(ks->f, ks->buf, __bufsize);		\
					if (ks->end == 0) { ks->is_eof = 1; break; }		\
					if (ks->end == -1) { ks->is_eof = 1; return -3; }	\
				} else break;											\
			}															\
			if (delimiter == KS_SEP_LINE) { \
				for (i = ks->begin; i < ks->end; ++i) \
					if (ks->buf[i] == '\n') break; \
			} else if (delimiter > KS_SEP_MAX) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (ks->buf[i] == delimiter) break;					\
			} else if (delimiter == KS_SEP_SPACE) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (isspace(ks->buf[i])) break;						\
			} else if (delimiter == KS_SEP_TAB) {						\
				for (i = ks->begin; i < ks->end; ++i)					\
					if (isspace(ks->buf[i]) && ks->buf[i] != ' ') break; \
			} else i = 0; /* never come to here! */						\
			if (str->m - str->l < (size_t)(i - ks->begin + 1)) {		\
				str->m = str->l + (i - ks->begin) + 1;					\
				kroundup32(str->m);										\
				str->s = (char*)realloc(str->s, str->m);				\
			}															\
			gotany = 1;													\
			memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin); \
			str->l = str->l + (i - ks->begin);							\
			ks->begin = i + 1;											\
			if (i < ks->end) {											\
				if (dret) *dret = ks->buf[i];							\
				break;													\
			}															\
		}																\
		if (!gotany && ks_eof(ks)) return -1;							\
		if (str->s == 0) {												\
			str->m = 1;													\
			str->s = (char*)calloc(1, 1);								\
		} else if (delimiter == KS_SEP_LINE && str->l > 1 && str->s[str->l-1] == '\r') --str->l; \
		str->s[str->l] = '\0';											\
		return str->l;													\
	} \
	int GLUE(ks_getuntil_, reader)(GLUE(kstream_t_, reader) *ks, int delimiter, kstring_t *str, int *dret) \
	{ return GLUE(ks_getuntil2_, reader)(ks, delimiter, str, dret, 0); }

#define KSTREAM_INIT(type_t, __read, __bufsize, reader) \
	__KS_TYPE(type_t, reader)							\
	__KS_BASIC(type_t, __bufsize, reader)				\
	__KS_GETC(__read, __bufsize, reader)				\
	__KS_GETUNTIL(__read, __bufsize, reader)

#define kseq_rewind(ks) ((ks)->last_char = (ks)->f->is_eof = (ks)->f->begin = (ks)->f->end = 0)

#define __KSEQ_BASIC(SCOPE, type_t, reader)										\
	SCOPE GLUE(kseq_t_, reader) *GLUE(kseq_init_, reader)(type_t fd)									\
	{																	\
		GLUE(kseq_t_, reader) *s = (GLUE(kseq_t_, reader)*)calloc(1, sizeof(GLUE(kseq_t_, reader)));					\
		s->f = GLUE(ks_init_, reader)(fd);												\
		return s;														\
	}																	\
	SCOPE void GLUE(kseq_destroy_, reader)(GLUE(kseq_t_, reader) *ks)									\
	{																	\
		if (!ks) return;												\
		free(ks->name.s); free(ks->comment.s); free(ks->seq.s);	free(ks->qual.s); \
		GLUE(ks_destroy_, reader)(ks->f);												\
		free(ks);														\
	}

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
   -3   error reading stream
 */
#define __KSEQ_READ(SCOPE, reader) \
	SCOPE int GLUE(kseq_read_, reader)(GLUE(kseq_t_, reader) *seq) \
	{ \
		int c,r; \
		GLUE(kstream_t_, reader) *ks = seq->f; \
		if (seq->last_char == 0) { /* then jump to the next header line */ \
			while ((c = GLUE(ks_getc_, reader)(ks)) >= 0 && c != '>' && c != '@'); \
			if (c < 0) return c; /* end of file or error*/ \
			seq->last_char = c; \
		} /* else: the first header char has been read in the previous call */ \
		seq->comment.l = seq->seq.l = seq->qual.l = 0; /* reset all members */ \
		if ((r=GLUE(ks_getuntil_, reader)(ks, 0, &seq->name, &c)) < 0) return r;  /* normal exit: EOF or error */ \
		if (c != '\n') GLUE(ks_getuntil_, reader)(ks, KS_SEP_LINE, &seq->comment, 0); /* read FASTA/Q comment */ \
		if (seq->seq.s == 0) { /* we can do this in the loop below, but that is slower */ \
			seq->seq.m = 256; \
			seq->seq.s = (char*)malloc(seq->seq.m); \
		} \
		while ((c = GLUE(ks_getc_, reader)(ks)) >= 0 && c != '>' && c != '+' && c != '@') { \
			if (c == '\n') continue; /* skip empty lines */ \
			seq->seq.s[seq->seq.l++] = c; /* this is safe: we always have enough space for 1 char */ \
			GLUE(ks_getuntil2_, reader)(ks, KS_SEP_LINE, &seq->seq, 0, 1); /* read the rest of the line */ \
		} \
		if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */	\
		if (seq->seq.l + 1 >= seq->seq.m) { /* seq->seq.s[seq->seq.l] below may be out of boundary */ \
			seq->seq.m = seq->seq.l + 2; \
			kroundup32(seq->seq.m); /* rounded to the next closest 2^k */ \
			seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m); \
		} \
		seq->seq.s[seq->seq.l] = 0;	/* null terminated string */ \
		if (c != '+') return seq->seq.l; /* FASTA */ \
		if (seq->qual.m < seq->seq.m) {	/* allocate memory for qual in case insufficient */ \
			seq->qual.m = seq->seq.m; \
			seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m); \
		} \
		while ((c = GLUE(ks_getc_, reader)(ks)) >= 0 && c != '\n'); /* skip the rest of '+' line */ \
		if (c == -1) return -2; /* error: no quality string */ \
		while ((c = GLUE(ks_getuntil2_, reader)(ks, KS_SEP_LINE, &seq->qual, 0, 1) >= 0 && seq->qual.l < seq->seq.l)); \
		if (c == -3) return -3; /* stream error */ \
		seq->last_char = 0;	/* we have not come to the next header line */ \
		if (seq->seq.l != seq->qual.l) return -2; /* error: qual string is of a different length */ \
		return seq->seq.l; \
	}

#define __KSEQ_TYPE(type_t, reader)						\
	typedef struct {							\
		kstring_t name, comment, seq, qual;		\
		int last_char;							\
		GLUE(kstream_t_, reader) *f;							\
	} GLUE(kseq_t_, reader);

#define KSEQ_INIT2(SCOPE, type_t, __read, reader)		\
	KSTREAM_INIT(type_t, __read, 16384, reader)			\
	__KSEQ_TYPE(type_t, reader)							\
	__KSEQ_BASIC(SCOPE, type_t, reader)					\
	__KSEQ_READ(SCOPE, reader)

#define KSEQ_INIT(type_t, __read) KSEQ_INIT2(static, type_t, __read, default_reader)
#define KSEQ_INIT_READER(type_t, __read, reader) KSEQ_INIT2(static, type_t, __read, reader)

