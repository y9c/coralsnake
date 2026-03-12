#include <Python.h>
#include <string.h>
#include <zlib.h>

static inline int fast_itoa(int val, char *buf) {
    if (val == 0) { *buf = '0'; return 1; }
    char tmp[12]; int len = 0;
    while (val > 0) { tmp[len++] = (val % 10) + '0'; val /= 10; }
    for (int i = 0; i < len / 2; i++) { char t = tmp[i]; tmp[i] = tmp[len-1-i]; tmp[len-1-i] = t; }
    memcpy(buf, tmp, len);
    return len;
}

static PyObject* score_and_tag(PyObject* self, PyObject* args) {
    const char *cigar_str, *seq, *ref; int is_o1;
    if (!PyArg_ParseTuple(args, "sssp", &cigar_str, &seq, &ref, &is_o1)) return NULL;
    Py_ssize_t q_len = strlen(seq), r_len = strlen(ref);
    size_t buf_sz = (size_t)q_len * 16 + (size_t)r_len + 128;
    char* md_buf = (char*)malloc(buf_sz);
    if (!md_buf) return PyErr_NoMemory();
    int yf=0, zf=0, yc=0, zc=0, ns=0, nc=0;
    Py_ssize_t r_idx=0, q_idx=0;
    int match_count=0; size_t md_pos=0;
    int matches=0, exp_conv=0, wr_conv=0, other_mm=0, indels=0;
    char b1 = is_o1 ? 'A' : 'T', b3 = is_o1 ? 'C' : 'G';
    long length = 0;
    for (int i = 0; cigar_str[i]; i++) {
        if (cigar_str[i] >= '0' && cigar_str[i] <= '9') length = length * 10 + (cigar_str[i] - '0');
        else {
            char op = cigar_str[i];
            if (op == 'M' || op == '=' || op == 'X') {
                for (long j = 0; j < length && r_idx < r_len && q_idx < q_len; j++) {
                    char rb = ref[r_idx], qb = seq[q_idx];
                    if (rb == qb) { matches++; match_count++; if (qb == b1) zf++; else if (qb == b3) zc++; }
                    else {
                        md_pos += fast_itoa(match_count, md_buf + md_pos);
                        md_buf[md_pos++] = rb; match_count = 0;
                        if (is_o1) { if ((rb == 'C' && qb == 'T') || (rb == 'A' && qb == 'G')) { exp_conv++; if (qb == 'G') yf++; else yc++; } else { wr_conv++; ns++; } }
                        else { if ((rb == 'G' && qb == 'A') || (rb == 'T' && qb == 'C')) { exp_conv++; if (qb == 'A') yf++; else yc++; } else { wr_conv++; ns++; } }
                    }
                    r_idx++; q_idx++;
                }
            } else if (op == 'I' || op == 'S') { q_idx += length; nc += (int)length; if (op == 'I') indels += (int)length; }
            else if (op == 'D') {
                md_pos += fast_itoa(match_count, md_buf + md_pos); md_buf[md_pos++] = '^';
                for (long j = 0; j < length && r_idx < r_len; j++) md_buf[md_pos++] = ref[r_idx++];
                match_count = 0; nc += (int)length; indels += (int)length;
            } else if (op == 'N') r_idx += length;
            length = 0;
        }
    }
    md_pos += fast_itoa(match_count, md_buf + md_pos);
    if (q_idx != q_len) { free(md_buf); return Py_BuildValue("(iiOiiiiii)", -999, 999, Py_None, 0, 0, 0, 0, 0, 0); }
    int score = matches + exp_conv - wr_conv - other_mm - indels;
    PyObject* res = Py_BuildValue("(iiNiiiiii)", score, wr_conv + other_mm, PyUnicode_FromStringAndSize(md_buf, (Py_ssize_t)md_pos), yf, zf, yc, zc, ns, nc);
    free(md_buf); return res;
}

static PyObject* reverse_complement(PyObject* self, PyObject* args) {
    const char* seq; if (!PyArg_ParseTuple(args, "s", &seq)) return NULL;
    Py_ssize_t len = strlen(seq); char* rc = (char*)malloc((size_t)len + 1);
    static char complement[256]; static int initialized = 0;
    if (!initialized) {
        for (int i = 0; i < 256; i++) complement[i] = (char)i;
        complement['A'] = 'T'; complement['C'] = 'G'; complement['G'] = 'C'; complement['T'] = 'A';
        complement['a'] = 't'; complement['c'] = 'g'; complement['g'] = 'c'; complement['t'] = 'a';
        initialized = 1;
    }
    for (Py_ssize_t i = 0; i < len; i++) rc[len - 1 - i] = complement[(unsigned char)seq[i]];
    rc[len] = '\0'; PyObject* result = PyUnicode_FromStringAndSize(rc, len); free(rc); return result;
}

static PyObject* batch_base_conversion(PyObject* self, PyObject* args) {
    PyObject* seq_list; const char *from_bases, *to_bases;
    if (!PyArg_ParseTuple(args, "O!ss", &PyList_Type, &seq_list, &from_bases, &to_bases)) return NULL;
    Py_ssize_t n = PyList_Size(seq_list); PyObject* result_list = PyList_New(n);
    char lookup[256]; memset(lookup, 0, 256);
    for (int i = 0; from_bases[i] && to_bases[i]; i++) {
        lookup[(unsigned char)from_bases[i]] = to_bases[i];
        if (from_bases[i] >= 'A' && from_bases[i] <= 'Z') lookup[(unsigned char)from_bases[i] + 32] = to_bases[i] + 32;
        else if (from_bases[i] >= 'a' && from_bases[i] <= 'z') lookup[(unsigned char)from_bases[i] - 32] = to_bases[i] - 32;
    }
    for (Py_ssize_t i = 0; i < n; i++) {
        Py_ssize_t len; const char* seq = PyUnicode_AsUTF8AndSize(PyList_GetItem(seq_list, i), &len);
        char* target = (char*)malloc((size_t)len + 1);
        for (Py_ssize_t j = 0; j < len; j++) { unsigned char c = (unsigned char)seq[j]; target[j] = lookup[c] ? lookup[c] : (char)c; }
        target[len] = '\0'; PyList_SetItem(result_list, i, PyUnicode_FromStringAndSize(target, len)); free(target);
    }
    return result_list;
}

static PyObject* convert_fasta_file(PyObject* self, PyObject* args) {
    const char *in_fn, *out_fa, *f_bases, *t_fa;
    if (!PyArg_ParseTuple(args, "ssss", &in_fn, &out_fa, &f_bases, &t_fa)) return NULL;
    gzFile in = gzopen(in_fn, "r"); if (!in) return PyErr_SetFromErrno(PyExc_IOError);
    FILE* out = fopen(out_fa, "w"); if (!out) { gzclose(in); return PyErr_SetFromErrno(PyExc_IOError); }
    char lookup[256]; memset(lookup, 0, 256);
    for (int i=0; f_bases[i] && t_fa[i]; i++) {
        lookup[(unsigned char)f_bases[i]] = t_fa[i];
        if (f_bases[i]>='A' && f_bases[i]<='Z') lookup[(unsigned char)f_bases[i]+32] = t_fa[i]+32;
        else if (f_bases[i]>='a' && f_bases[i]<='z') lookup[(unsigned char)f_bases[i]-32] = t_fa[i]-32;
    }
    char buf[65536];
    while (gzgets(in, buf, sizeof(buf))) {
        if (buf[0] == '>') fputs(buf, out);
        else { for (int i=0; buf[i]; i++) { unsigned char c = (unsigned char)buf[i]; if (lookup[c]) buf[i] = lookup[c]; } fputs(buf, out); }
    }
    gzclose(in); fclose(out); Py_RETURN_NONE;
}

static PyMethodDef SeqOpsMethods[] = {
    {"score_and_tag", score_and_tag, METH_VARARGS, ""},
    {"reverse_complement", reverse_complement, METH_VARARGS, ""},
    {"batch_base_conversion", batch_base_conversion, METH_VARARGS, ""},
    {"convert_fasta_file", convert_fasta_file, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef seqopsmodule = { PyModuleDef_HEAD_INIT, "seqops", "", -1, SeqOpsMethods };
PyMODINIT_FUNC PyInit_seqops(void) { return PyModule_Create(&seqopsmodule); }
