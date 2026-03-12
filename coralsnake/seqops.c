#include <Python.h>
#include <string.h>
#include <zlib.h>

static PyObject* score_and_tag(PyObject* self, PyObject* args) {
    const char *cigar_str, *seq, *ref;
    int is_o1;
    if (!PyArg_ParseTuple(args, "sssp", &cigar_str, &seq, &ref, &is_o1)) return NULL;
    Py_ssize_t q_len = strlen(seq), r_len = strlen(ref);
    size_t buf_sz = (size_t)q_len * 16 + (size_t)r_len + 128; // Ensure space for massive D ops
    char* md_buf = (char*)malloc(buf_sz);
    if (!md_buf) return PyErr_NoMemory();
    int yf=0, zf=0, yc=0, zc=0, ns=0, nc=0;
    Py_ssize_t r_idx=0, q_idx=0;
    int match_count=0;
    size_t md_pos=0;
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
                        int written = snprintf(md_buf + md_pos, buf_sz - md_pos, "%d%c", match_count, rb);
                        if (written > 0 && (size_t)written < buf_sz - md_pos) md_pos += written;
                        match_count = 0;
                        if (is_o1) {
                            if ((rb == 'C' && qb == 'T') || (rb == 'A' && qb == 'G')) { exp_conv++; if (qb == 'G') yf++; else yc++; }
                            else { wr_conv++; ns++; }
                        } else {
                            if ((rb == 'G' && qb == 'A') || (rb == 'T' && qb == 'C')) { exp_conv++; if (qb == 'A') yf++; else yc++; }
                            else { wr_conv++; ns++; }
                        }
                    }
                    r_idx++; q_idx++;
                }
            } else if (op == 'I' || op == 'S') { q_idx += length; nc += (int)length; if (op == 'I') indels += (int)length; }
            else if (op == 'D') {
                int written = snprintf(md_buf + md_pos, buf_sz - md_pos, "%d^", match_count);
                if (written > 0 && (size_t)written < buf_sz - md_pos) md_pos += written;
                for (long j = 0; j < length && r_idx < r_len; j++) {
                    if (md_pos < buf_sz - 1) md_buf[md_pos++] = ref[r_idx++];
                    else r_idx++;
                }
                match_count = 0; nc += (int)length; indels += (int)length;
            } else if (op == 'N') r_idx += length;
            length = 0;
        }
    }
    int written = snprintf(md_buf + md_pos, buf_sz - md_pos, "%d", match_count);
    if (written > 0 && (size_t)written < buf_sz - md_pos) md_pos += written;
    
    if (q_idx != q_len) { free(md_buf); return Py_BuildValue("(iiOiiiiii)", -999, 999, Py_None, 0, 0, 0, 0, 0, 0); }
    int score = matches + exp_conv - wr_conv - other_mm - indels;
    PyObject* res = Py_BuildValue("(iiNiiiiii)", score, wr_conv + other_mm, PyUnicode_FromStringAndSize(md_buf, (Py_ssize_t)md_pos), yf, zf, yc, zc, ns, nc);
    free(md_buf); return res;
}

static PyObject* batch_score_and_tag(PyObject* self, PyObject* args) {
    PyObject *hit_list, *q_orig_list, *ref_list;
    int is_o1;
    if (!PyArg_ParseTuple(args, "OOOp", &hit_list, &q_orig_list, &ref_list, &is_o1)) return NULL;
    Py_ssize_t n = PyList_Size(hit_list);
    PyObject* result_list = PyList_New(n);
    for (Py_ssize_t i = 0; i < n; i++) {
        PyObject* hit = PyList_GetItem(hit_list, i);
        PyObject* q_item = PyList_GetItem(q_orig_list, i);
        PyObject* r_item = PyList_GetItem(ref_list, i);
        
        if (hit == Py_None || q_item == Py_None || r_item == Py_None) {
            Py_INCREF(Py_None); PyList_SetItem(result_list, i, Py_None); continue;
        }
        
        PyObject* cigar_obj = PyTuple_GetItem(hit, 7);
        const char* cigar = PyUnicode_AsUTF8(cigar_obj);
        const char* q_orig = PyUnicode_AsUTF8(q_item);
        const char* ref = PyUnicode_AsUTF8(r_item);
        
        if (!q_orig || !ref || !cigar) {
            PyErr_Clear(); // IMPORTANT: Clear error from AsUTF8 if any
            Py_INCREF(Py_None); PyList_SetItem(result_list, i, Py_None); continue;
        }
        
        PyObject* score_args = Py_BuildValue("(sssp)", cigar, q_orig, ref, is_o1);
        PyObject* scored = score_and_tag(NULL, score_args);
        Py_DECREF(score_args);
        PyList_SetItem(result_list, i, scored);
    }
    return result_list;
}

static PyObject* reverse_complement(PyObject* self, PyObject* args) {
    const char* seq;
    if (!PyArg_ParseTuple(args, "s", &seq)) return NULL;
    Py_ssize_t len = strlen(seq);
    char* rc = (char*)malloc((size_t)len + 1);
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
    {"batch_score_and_tag", batch_score_and_tag, METH_VARARGS, ""},
    {"reverse_complement", reverse_complement, METH_VARARGS, ""},
    {"batch_base_conversion", batch_base_conversion, METH_VARARGS, ""},
    {"convert_fasta_file", convert_fasta_file, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef seqopsmodule = { PyModuleDef_HEAD_INIT, "seqops", "", -1, SeqOpsMethods };
PyMODINIT_FUNC PyInit_seqops(void) { return PyModule_Create(&seqopsmodule); }
