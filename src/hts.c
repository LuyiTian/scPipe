/*  hts.c -- format-neutral I/O, indexing, and iterator API functions.

    Copyright (C) 2008, 2009, 2012-2015 Genome Research Ltd.
    Copyright (C) 2012, 2013 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */



#include <zlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/stat.h>
#include "bgzf.h"
#include "hts.h"
#include "hfile.h"
#include "hts_internal.h"

#include "kstring.h"
#include "kseq.h"
#define KS_BGZF 1
#if KS_BGZF
    // bgzf now supports gzip-compressed files, the gzFile branch can be removed
    KSTREAM_INIT2(, BGZF*, bgzf_read, 65536)
#endif

#include "khash.h"
KHASH_INIT2(s2i,, kh_cstr_t, int64_t, 1, kh_str_hash_func, kh_str_hash_equal)

int hts_verbose = 3;


const unsigned char seq_nt16_table[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
     1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";

const int seq_nt16_int[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

/**********************
 *** Basic file I/O ***
 **********************/

static enum htsFormatCategory format_category(enum htsExactFormat fmt)
{
    switch (fmt) {
    case bam:
    case sam:

    case vcf:
    case bcf:
        return variant_data;

    case bai:
    case crai:
    case csi:
    case gzi:
    case tbi:
        return index_file;

    case bed:
        return region_list;

    case unknown_format:
    case binary_format:
    case text_format:
    case format_maximum:
        break;
    }

    return unknown_category;
}

// Decompress up to ten or so bytes by peeking at the file, which must be
// positioned at the start of a GZIP block.
static size_t decompress_peek(hFILE *fp, unsigned char *dest, size_t destsize)
{
    // Typically at most a couple of hundred bytes of input are required
    // to get a few bytes of output from inflate(), so hopefully this buffer
    // size suffices in general.
    unsigned char buffer[512];
    z_stream zs;
    ssize_t npeek = hpeek(fp, buffer, sizeof buffer);

    if (npeek < 0) return 0;

    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = buffer;
    zs.avail_in = npeek;
    zs.next_out = dest;
    zs.avail_out = destsize;
    if (inflateInit2(&zs, 31) != Z_OK) return 0;

    while (zs.total_out < destsize)
        if (inflate(&zs, Z_SYNC_FLUSH) != Z_OK) break;

    destsize = zs.total_out;
    inflateEnd(&zs);

    return destsize;
}

// Parse "x.y" text, taking care because the string is not NUL-terminated
// and filling in major/minor only when the digits are followed by a delimiter,
// so we don't misread "1.10" as "1.1" due to reaching the end of the buffer.
static void
parse_version(htsFormat *fmt, const unsigned char *u, const unsigned char *ulim)
{
    const char *str  = (const char *) u;
    const char *slim = (const char *) ulim;
    const char *s;

    fmt->version.major = fmt->version.minor = -1;

    for (s = str; s < slim; s++) if (!isdigit(*s)) break;
    if (s < slim) {
        fmt->version.major = atoi(str);
        if (*s == '.') {
            str = &s[1];
            for (s = str; s < slim; s++) if (!isdigit(*s)) break;
            if (s < slim)
                fmt->version.minor = atoi(str);
        }
        else
            fmt->version.minor = 0;
    }
}

int hts_detect_format(hFILE *hfile, htsFormat *fmt)
{
    unsigned char s[21];
    ssize_t len = hpeek(hfile, s, 18);
    if (len < 0) return -1;

    if (len >= 2 && s[0] == 0x1f && s[1] == 0x8b) {
        // The stream is either gzip-compressed or BGZF-compressed.
        // Determine which, and decompress the first few bytes.
        fmt->compression = (len >= 18 && (s[3] & 4) &&
                            memcmp(&s[12], "BC\2\0", 4) == 0)? bgzf : gzip;
        len = decompress_peek(hfile, s, sizeof s);
    }
    else {
        fmt->compression = no_compression;
        len = hpeek(hfile, s, sizeof s);
    }
    if (len < 0) return -1;

    fmt->compression_level = -1;
    fmt->specific = NULL;


    if (len >= 4 && s[3] <= '\4') {
        if (memcmp(s, "BAM\1", 4) == 0) {
            fmt->category = sequence_data;
            fmt->format = bam;
            // TODO Decompress enough to pick version from @HD-VN header
            fmt->version.major = 1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "BAI\1", 4) == 0) {
            fmt->category = index_file;
            fmt->format = bai;
            fmt->version.major = -1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "BCF\4", 4) == 0) {
            fmt->category = variant_data;
            fmt->format = bcf;
            fmt->version.major = 1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "BCF\2", 4) == 0) {
            fmt->category = variant_data;
            fmt->format = bcf;
            fmt->version.major = s[3];
            fmt->version.minor = (len >= 5 && s[4] <= 2)? s[4] : 0;
            return 0;
        }
        else if (memcmp(s, "CSI\1", 4) == 0) {
            fmt->category = index_file;
            fmt->format = csi;
            fmt->version.major = 1, fmt->version.minor = -1;
            return 0;
        }
        else if (memcmp(s, "TBI\1", 4) == 0) {
            fmt->category = index_file;
            fmt->format = tbi;
            fmt->version.major = -1, fmt->version.minor = -1;
            return 0;
        }
    }
    else if (len >= 16 && memcmp(s, "##fileformat=VCF", 16) == 0) {
        fmt->category = variant_data;
        fmt->format = vcf;
        if (len >= 21 && s[16] == 'v')
            parse_version(fmt, &s[17], &s[len]);
        else
            fmt->version.major = fmt->version.minor = -1;
        return 0;
    }
    else if (len >= 4 && s[0] == '@' &&
             (memcmp(s, "@HD\t", 4) == 0 || memcmp(s, "@SQ\t", 4) == 0 ||
              memcmp(s, "@RG\t", 4) == 0 || memcmp(s, "@PG\t", 4) == 0)) {
        fmt->category = sequence_data;
        fmt->format = sam;
        // @HD-VN is not guaranteed to be the first tag, but then @HD is
        // not guaranteed to be present at all...
        if (len >= 9 && memcmp(s, "@HD\tVN:", 7) == 0)
            parse_version(fmt, &s[7], &s[len]);
        else
            fmt->version.major = 1, fmt->version.minor = -1;
        return 0;
    }
    else {
        // Various possibilities for tab-delimited text:
        // .crai   (gzipped tab-delimited six columns: seqid 5*number)
        // .bed    ([3..12] tab-delimited columns)
        // .bedpe  (>= 10 tab-delimited columns)
        // .sam    (tab-delimited >= 11 columns: seqid number seqid...)
        // FIXME For now, assume it's SAM
        fmt->category = sequence_data;
        fmt->format = sam;
        fmt->version.major = 1, fmt->version.minor = -1;
        return 0;
    }

    fmt->category = unknown_category;
    fmt->format = unknown_format;
    fmt->version.major = fmt->version.minor = -1;
    fmt->compression = no_compression;
    return 0;
}

char *hts_format_description(const htsFormat *format)
{
    kstring_t str = { 0, 0, NULL };

    switch (format->format) {
    case sam:   kputs("SAM", &str); break;
    case bam:   kputs("BAM", &str); break;
    case vcf:   kputs("VCF", &str); break;
    case bcf:
        if (format->version.major == 1) kputs("Legacy BCF", &str);
        else kputs("BCF", &str);
        break;
    case bai:   kputs("BAI", &str); break;
    case crai:  kputs("CRAI", &str); break;
    case csi:   kputs("CSI", &str); break;
    case tbi:   kputs("Tabix", &str); break;
    default:    kputs("unknown", &str); break;
    }

    if (format->version.major >= 0) {
        kputs(" version ", &str);
        kputw(format->version.major, &str);
        if (format->version.minor >= 0) {
            kputc('.', &str);
            kputw(format->version.minor, &str);
        }
    }

    switch (format->compression) {
    case custom: kputs(" compressed", &str); break;
    case gzip:   kputs(" gzip-compressed", &str); break;
    case bgzf:
        switch (format->format) {
        case bam:
        case bcf:
        case csi:
        case tbi:
            // These are by definition BGZF, so just use the generic term
            kputs(" compressed", &str);
            break;
        default:
            kputs(" BGZF-compressed", &str);
            break;
        }
        break;
    default: break;
    }

    switch (format->category) {
    case sequence_data: kputs(" sequence", &str); break;
    case variant_data:  kputs(" variant calling", &str); break;
    case index_file:    kputs(" index", &str); break;
    case region_list:   kputs(" genomic region", &str); break;
    default: break;
    }

    if (format->compression == no_compression)
        switch (format->format) {
        case sam:
        case crai:
        case vcf:
        case bed:
            kputs(" text", &str);
            break;

        default:
            kputs(" data", &str);
            break;
        }
    else
        kputs(" data", &str);

    return ks_release(&str);
}

htsFile *hts_open_format(const char *fn, const char *mode, const htsFormat *fmt)
{
    char smode[102], *cp, *cp2, *mode_c;
    htsFile *fp = NULL;
    hFILE *hfile;
    char fmt_code = '\0';

    strncpy(smode, mode, 100);
    smode[100]=0;
    if ((cp = strchr(smode, ',')))
        *cp = '\0';

    // Migrate format code (b or c) to the end of the smode buffer.
    for (cp2 = cp = smode; *cp; cp++) {
        if (*cp == 'b')
            fmt_code = 'b';
        else if (*cp == 'c')
            fmt_code = 'c';
        else
            *cp2++ = *cp;
    }
    mode_c = cp2;
    *cp2++ = fmt_code;
    *cp2++ = 0;
    *cp2++ = 0;

    // Set or reset the format code if opts->format is used
    if (fmt && fmt->format != unknown_format)
        *mode_c = "\0g\0\0b\0c\0\0b\0g\0\0"[fmt->format];

    hfile = hopen(fn, smode);
    if (hfile == NULL) goto error;

    fp = hts_hopen(hfile, fn, smode);
    if (fp == NULL) goto error;


    return fp;

error:
    if (hts_verbose >= 2)
        fprintf(stderr, "[E::%s] fail to open file '%s'\n", __func__, fn);

    if (hfile)
        hclose_abruptly(hfile);

    return NULL;
}

htsFile *hts_open(const char *fn, const char *mode) {
    return hts_open_format(fn, mode, NULL);
}

/*
 * Parses arg and appends it to the option list.
 *
 * Returns 0 on success;
 *        -1 on failure.
 */
int hts_opt_add(hts_opt **opts, const char *c_arg) {
    hts_opt *o, *t;
    char *val;

    if (!c_arg)
        return -1;

    if (!(o =  malloc(sizeof(*o))))
        return -1;

    if (!(o->arg = strdup(c_arg))) {
        free(o);
        return -1;
    }

    if (!(val = strchr(o->arg, '=')))
        val = "1"; // assume boolean
    else
        *val++ = '\0';


    o->next = NULL;

    // Append; assumes small list.
    if (*opts) {
        t = *opts;
        while (t->next)
            t = t->next;
        t->next = o;
    } else {
        *opts = o;
    }

    return 0;
}

/*
 * Applies an hts_opt option list to a given htsFile.
 *
 * Returns 0 on success
 *        -1 on failure
 */
int hts_opt_apply(htsFile *fp, hts_opt *opts) {
    hts_opt *last = NULL;

    for (; opts;  opts = (last=opts)->next)
        if (hts_set_opt(fp,  opts->opt,  opts->val) != 0)
            return -1;

    return 0;
}

/*
 * Frees an hts_opt list.
 */
void hts_opt_free(hts_opt *opts) {
    hts_opt *last = NULL;
    while (opts) {
        opts = (last=opts)->next;
        free(last->arg);
        free(last);
    }
}


/*
 * Tokenise options as (key(=value)?,)*(key(=value)?)?
 * NB: No provision for ',' appearing in the value!
 * Add backslashing rules?
 *
 * This could be used as part of a general command line option parser or
 * as a string concatenated onto the file open mode.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
int hts_parse_opt_list(htsFormat *fmt, const char *str) {
    while (str && *str) {
        const char *str_start;
        int len;
        char arg[8001];

        while (*str && *str == ',')
            str++;

        for (str_start = str; *str && *str != ','; str++);
        len = str - str_start;

        // Produce a nul terminated copy of the option
        strncpy(arg, str_start, len < 8000 ? len : 8000);
        arg[len < 8000 ? len : 8000] = '\0';

        if (hts_opt_add((hts_opt **)&fmt->specific, arg) != 0)
            return -1;

        if (*str)
            str++;
    }

    return 0;
}

/*
 * Accepts a string file format (sam, bam, vcf, bam) optionally
 * followed by a comma separated list of key=value options and splits
 * these up into the fields of htsFormat struct.
 *
 * format is assumed to be already initialised, either to blank
 * "unknown" values or via previous hts_opt_add calls.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
int hts_parse_format(htsFormat *format, const char *str) {
    const char *cp;

    if (!(cp = strchr(str, ',')))
        cp = str+strlen(str);

    format->version.minor = 0; // unknown
    format->version.major = 0; // unknown

    if (strncmp(str, "sam", cp-str) == 0) {
        format->category          = sequence_data;
        format->format            = sam;
        format->compression       = no_compression;;
        format->compression_level = 0;
    } else if (strncmp(str, "bam", cp-str) == 0) {
        format->category          = sequence_data;
        format->format            = bam;
        format->compression       = bgzf;
        format->compression_level = -1;
    }  else if (strncmp(str, "vcf", cp-str) == 0) {
        format->category          = variant_data;
        format->format            = vcf;
        format->compression       = no_compression;;
        format->compression_level = 0;
    } else if (strncmp(str, "bcf", cp-str) == 0) {
        format->category          = variant_data;
        format->format            = bcf;
        format->compression       = bgzf;
        format->compression_level = -1;
    } else {
        return -1;
    }

    return hts_parse_opt_list(format, cp);
}


/*
 * Tokenise options as (key(=value)?,)*(key(=value)?)?
 * NB: No provision for ',' appearing in the value!
 * Add backslashing rules?
 *
 * This could be used as part of a general command line option parser or
 * as a string concatenated onto the file open mode.
 *
 * Returns 0 on success
 *        -1 on failure.
 */
static int hts_process_opts(htsFile *fp, const char *opts) {
    htsFormat fmt;

    fmt.specific = NULL;
    if (hts_parse_opt_list(&fmt, opts) != 0)
        return -1;

    if (hts_opt_apply(fp, fmt.specific) != 0) {
        hts_opt_free(fmt.specific);
        return -1;
    }

    hts_opt_free(fmt.specific);

    return 0;
}


htsFile *hts_hopen(struct hFILE *hfile, const char *fn, const char *mode)
{
    htsFile *fp = (htsFile*)calloc(1, sizeof(htsFile));
    char simple_mode[101], *cp, *opts;
    simple_mode[100] = '\0';

    if (fp == NULL) goto error;

    fp->fn = strdup(fn);
    fp->is_be = ed_is_big();

    // Split mode into simple_mode,opts strings
    if ((cp = strchr(mode, ','))) {
        strncpy(simple_mode, mode, cp-mode <= 100 ? cp-mode : 100);
        simple_mode[cp-mode] = '\0';
        opts = cp+1;
    } else {
        strncpy(simple_mode, mode, 100);
        opts = NULL;
    }

    if (strchr(simple_mode, 'r')) {
        if (hts_detect_format(hfile, &fp->format) < 0) goto error;
    }
    else if (strchr(simple_mode, 'w') || strchr(simple_mode, 'a')) {
        htsFormat *fmt = &fp->format;
        fp->is_write = 1;

        if (strchr(simple_mode, 'b')) fmt->format = binary_format;
        else fmt->format = text_format;

        if (strchr(simple_mode, 'z')) fmt->compression = bgzf;
        else if (strchr(simple_mode, 'g')) fmt->compression = gzip;
        else if (strchr(simple_mode, 'u')) fmt->compression = no_compression;
        else {
            // No compression mode specified, set to the default for the format
            switch (fmt->format) {
            case binary_format: fmt->compression = bgzf; break;
            case text_format: fmt->compression = no_compression; break;
            default: abort();
            }
        }

        // Fill in category (if determinable; e.g. 'b' could be BAM or BCF)
        fmt->category = format_category(fmt->format);

        fmt->version.major = fmt->version.minor = -1;
        fmt->compression_level = -1;
        fmt->specific = NULL;
    }
    else goto error;

    switch (fp->format.format) {
    case binary_format:
    case bam:
    case bcf:
        fp->fp.bgzf = bgzf_hopen(hfile, simple_mode);
        if (fp->fp.bgzf == NULL) goto error;
        fp->is_bin = 1;
        break;


    case text_format:
    case sam:
    case vcf:
        if (!fp->is_write) {
        #if KS_BGZF
            BGZF *gzfp = bgzf_hopen(hfile, simple_mode);
        #else
            // TODO Implement gzip hFILE adaptor
            hclose(hfile); // This won't work, especially for stdin
            gzFile gzfp = strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
        #endif
            if (gzfp) fp->fp.voidp = ks_init(gzfp);
            else goto error;
        }
        else if (fp->format.compression != no_compression) {
            fp->fp.bgzf = bgzf_hopen(hfile, simple_mode);
            if (fp->fp.bgzf == NULL) goto error;
        }
        else
            fp->fp.hfile = hfile;
        break;

    default:
        goto error;
    }

    if (opts)
        hts_process_opts(fp, opts);

    return fp;

error:
    if (hts_verbose >= 2)
        fprintf(stderr, "[E::%s] fail to open file '%s'\n", __func__, fn);

    if (fp) {
        free(fp->fn);
        free(fp->fn_aux);
        free(fp);
    }
    return NULL;
}

int hts_close(htsFile *fp)
{
    int ret, save;

    switch (fp->format.format) {
    case binary_format:
    case bam:
    case bcf:
        ret = bgzf_close(fp->fp.bgzf);
        break;


    case text_format:
    case sam:
    case vcf:
        if (!fp->is_write) {
        #if KS_BGZF
            BGZF *gzfp = ((kstream_t*)fp->fp.voidp)->f;
            ret = bgzf_close(gzfp);
        #else
            gzFile gzfp = ((kstream_t*)fp->fp.voidp)->f;
            ret = gzclose(gzfp);
        #endif
            ks_destroy((kstream_t*)fp->fp.voidp);
        }
        else if (fp->format.compression != no_compression)
            ret = bgzf_close(fp->fp.bgzf);
        else
            ret = hclose(fp->fp.hfile);
        break;

    default:
        ret = -1;
        break;
    }

    save = errno;
    free(fp->fn);
    free(fp->fn_aux);
    free(fp->line.s);
    free(fp);
    errno = save;
    return ret;
}

const htsFormat *hts_get_format(htsFile *fp)
{
    return fp? &fp->format : NULL;
}

const char *hts_format_file_extension(const htsFormat *format) {
    if (!format)
        return "?";

    switch (format->format) {
    case sam:  return "sam";
    case bam:  return "bam";
    case bai:  return "bai";
    case crai: return "crai";
    case vcf:  return "vcf";
    case bcf:  return "bcf";
    case csi:  return "csi";
    case gzi:  return "gzi";
    case tbi:  return "tbi";
    case bed:  return "bed";
    default:   return "?";
    }
}

int hts_set_opt(htsFile *fp, enum hts_fmt_option opt, ...) {
    va_list args;

    if (opt == HTS_OPT_NTHREADS) {
        va_start(args, opt);
        int nthreads = 4;
        va_end(args);
        return hts_set_threads(fp, nthreads);
    }
    else
    {
        return 1;
    }

}

int hts_set_threads(htsFile *fp, int n)
{
    if (fp->format.compression == bgzf) {
        return bgzf_mt(fp->fp.bgzf, n, 256);
    }
    else return 0;
}

int hts_set_fai_filename(htsFile *fp, const char *fn_aux)
{
    free(fp->fn_aux);
    if (fn_aux) {
        fp->fn_aux = strdup(fn_aux);
        if (fp->fn_aux == NULL) return -1;
    }
    else fp->fn_aux = NULL;


    return 0;
}

// For VCF/BCF backward sweeper. Not exposing these functions because their
// future is uncertain. Things will probably have to change with hFILE...
BGZF *hts_get_bgzfp(htsFile *fp)
{
    if ( fp->is_bin )
        return fp->fp.bgzf;
    else
        return ((kstream_t*)fp->fp.voidp)->f;
}
int hts_useek(htsFile *fp, long uoffset, int where)
{
    if ( fp->is_bin )
        return bgzf_useek(fp->fp.bgzf, uoffset, where);
    else
    {
        ks_rewind((kstream_t*)fp->fp.voidp);
        ((kstream_t*)fp->fp.voidp)->seek_pos = uoffset;
        return bgzf_useek(((kstream_t*)fp->fp.voidp)->f, uoffset, where);
    }
}
long hts_utell(htsFile *fp)
{
    if ( fp->is_bin )
        return bgzf_utell(fp->fp.bgzf);
    else
        return ((kstream_t*)fp->fp.voidp)->seek_pos;
}

int hts_getline(htsFile *fp, int delimiter, kstring_t *str)
{
    int ret, dret;
    ret = ks_getuntil((kstream_t*)fp->fp.voidp, delimiter, str, &dret);
    ++fp->lineno;
    return ret;
}

char **hts_readlist(const char *string, int is_file, int *_n)
{
    int m = 0, n = 0, dret;
    char **s = 0;
    if ( is_file )
    {
#if KS_BGZF
        BGZF *fp = bgzf_open(string, "r");
#else
        gzFile fp = gzopen(string, "r");
#endif
        if ( !fp ) return NULL;

        kstream_t *ks;
        kstring_t str;
        str.s = 0; str.l = str.m = 0;
        ks = ks_init(fp);
        while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0)
        {
            if (str.l == 0) continue;
            n++;
            hts_expand(char*,n,m,s);
            s[n-1] = strdup(str.s);
        }
        ks_destroy(ks);
#if KS_BGZF
        bgzf_close(fp);
#else
        gzclose(fp);
#endif
        free(str.s);
    }
    else
    {
        const char *q = string, *p = string;
        while ( 1 )
        {
            if (*p == ',' || *p == 0)
            {
                n++;
                hts_expand(char*,n,m,s);
                s[n-1] = (char*)calloc(p - q + 1, 1);
                strncpy(s[n-1], q, p - q);
                q = p + 1;
            }
            if ( !*p ) break;
            p++;
        }
    }
    s = (char**)realloc(s, n * sizeof(char*));
    *_n = n;
    return s;
}

char **hts_readlines(const char *fn, int *_n)
{
    int m = 0, n = 0, dret;
    char **s = 0;
#if KS_BGZF
    BGZF *fp = bgzf_open(fn, "r");
#else
    gzFile fp = gzopen(fn, "r");
#endif
    if ( fp ) { // read from file
        kstream_t *ks;
        kstring_t str;
        str.s = 0; str.l = str.m = 0;
        ks = ks_init(fp);
        while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
            if (str.l == 0) continue;
            if (m == n) {
                m = m? m<<1 : 16;
                s = (char**)realloc(s, m * sizeof(char*));
            }
            s[n++] = strdup(str.s);
        }
        ks_destroy(ks);
        #if KS_BGZF
            bgzf_close(fp);
        #else
            gzclose(fp);
        #endif
        s = (char**)realloc(s, n * sizeof(char*));
        free(str.s);
    } else if (*fn == ':') { // read from string
        const char *q, *p;
        for (q = p = fn + 1;; ++p)
            if (*p == ',' || *p == 0) {
                if (m == n) {
                    m = m? m<<1 : 16;
                    s = (char**)realloc(s, m * sizeof(char*));
                }
                s[n] = (char*)calloc(p - q + 1, 1);
                strncpy(s[n++], q, p - q);
                q = p + 1;
                if (*p == 0) break;
            }
    } else return 0;
    s = (char**)realloc(s, n * sizeof(char*));
    *_n = n;
    return s;
}

// DEPRECATED: To be removed in a future HTSlib release
int hts_file_type(const char *fname)
{
    int len = strlen(fname);
    if ( !strcasecmp(".vcf.gz",fname+len-7) ) return FT_VCF_GZ;
    if ( !strcasecmp(".vcf",fname+len-4) ) return FT_VCF;
    if ( !strcasecmp(".bcf",fname+len-4) ) return FT_BCF_GZ;
    if ( !strcmp("-",fname) ) return FT_STDIN;

    hFILE *f = hopen(fname, "r");
    if (f == NULL) return 0;

    htsFormat fmt;
    if (hts_detect_format(f, &fmt) < 0) { hclose_abruptly(f); return 0; }
    if (hclose(f) < 0) return 0;

    switch (fmt.format) {
    case vcf: return (fmt.compression == no_compression)? FT_VCF : FT_VCF_GZ;
    case bcf: return (fmt.compression == no_compression)? FT_BCF : FT_BCF_GZ;
    default:  return 0;
    }
}


