#include <Python.h>
#include <string.h>

// Fast base conversion function
static PyObject* fast_base_conversion(PyObject* self, PyObject* args) {
    const char* seq;
    const char* from_bases;
    const char* to_bases;
    Py_ssize_t seq_len;
    
    if (!PyArg_ParseTuple(args, "sss", &seq, &from_bases, &to_bases)) {
        return NULL;
    }
    
    seq_len = strlen(seq);
    PyObject* result = PyUnicode_New(seq_len, 127); // ASCII max
    if (!result) {
        return NULL;
    }
    
    char* result_str = PyUnicode_AsUTF8(result);
    if (!result_str) {
        Py_DECREF(result);
        return NULL;
    }
    
    // Create lookup table for conversion
    char lookup[256];
    memset(lookup, 0, 256);
    
    // Build lookup table
    for (int i = 0; from_bases[i] && to_bases[i]; i++) {
        lookup[(unsigned char)from_bases[i]] = to_bases[i];
        // Also handle lowercase
        if (from_bases[i] >= 'A' && from_bases[i] <= 'Z') {
            lookup[(unsigned char)from_bases[i] + 32] = to_bases[i] + 32;
        } else if (from_bases[i] >= 'a' && from_bases[i] <= 'z') {
            lookup[(unsigned char)from_bases[i] - 32] = to_bases[i] - 32;
        }
    }
    
    // Convert sequence
    for (Py_ssize_t i = 0; i < seq_len; i++) {
        unsigned char c = (unsigned char)seq[i];
        result_str[i] = lookup[c] ? lookup[c] : c;
    }
    
    return result;
}

// Fast reverse complement function
static PyObject* fast_reverse_complement(PyObject* self, PyObject* args) {
    const char* seq;
    Py_ssize_t seq_len;
    
    if (!PyArg_ParseTuple(args, "s", &seq)) {
        return NULL;
    }
    
    seq_len = strlen(seq);
    PyObject* result = PyUnicode_New(seq_len, 127); // ASCII max
    if (!result) {
        return NULL;
    }
    
    char* result_str = PyUnicode_AsUTF8(result);
    if (!result_str) {
        Py_DECREF(result);
        return NULL;
    }
    
    // Complement lookup table
    char complement[256];
    memset(complement, 0, 256);
    complement['A'] = 'T';
    complement['T'] = 'A';
    complement['G'] = 'C';
    complement['C'] = 'G';
    complement['N'] = 'N';
    complement['a'] = 't';
    complement['t'] = 'a';
    complement['g'] = 'c';
    complement['c'] = 'g';
    complement['n'] = 'n';
    
    // Reverse and complement
    for (Py_ssize_t i = 0; i < seq_len; i++) {
        unsigned char c = (unsigned char)seq[seq_len - 1 - i];
        result_str[i] = complement[c] ? complement[c] : c;
    }
    
    return result;
}

// Fast MD tag and conversion statistics calculation
static PyObject* fast_cal_md_and_tag(PyObject* self, PyObject* args) {
    PyObject* cigar_list;
    const char* seq;
    const char* ref;
    int fwd;
    
    if (!PyArg_ParseTuple(args, "Ossp", &cigar_list, &seq, &ref, &fwd)) {
        return NULL;
    }
    
    Py_ssize_t cigar_len = PyList_Size(cigar_list);
    Py_ssize_t seq_len = strlen(seq);
    Py_ssize_t ref_len = strlen(ref);
    
    // Allocate MD tag buffer (worst case: every base is a mismatch)
    char* md_buffer = (char*)malloc(seq_len * 10);
    if (!md_buffer) {
        return PyErr_NoMemory();
    }
    
    int yf = 0, zf = 0, yc = 0, zc = 0, ns = 0, nc = 0;
    int ref_index = 0, query_index = 0, match_count = 0;
    int md_pos = 0;
    
    // Set bases based on direction
    char b1, b2, b3, b4;
    if (fwd) {
        b1 = 'A'; b2 = 'G'; b3 = 'C'; b4 = 'T';
    } else {
        b1 = 'T'; b2 = 'C'; b3 = 'G'; b4 = 'A';
    }
    
    // Process CIGAR operations
    for (Py_ssize_t i = 0; i < cigar_len; i++) {
        PyObject* cigar_op = PyList_GetItem(cigar_list, i);
        long length, operation;
        
        // Handle both list and tuple formats
        if (PyList_Check(cigar_op) && PyList_Size(cigar_op) == 2) {
            length = PyLong_AsLong(PyList_GetItem(cigar_op, 0));
            operation = PyLong_AsLong(PyList_GetItem(cigar_op, 1));
        } else if (PyTuple_Check(cigar_op) && PyTuple_Size(cigar_op) == 2) {
            length = PyLong_AsLong(PyTuple_GetItem(cigar_op, 0));
            operation = PyLong_AsLong(PyTuple_GetItem(cigar_op, 1));
        } else {
            free(md_buffer);
            PyErr_SetString(PyExc_ValueError, "CIGAR must be list of [length, operation] or (length, operation)");
            return NULL;
        }
        
        if (operation == 0) {  // Match or Mismatch
            for (long j = 0; j < length; j++) {
                if (ref_index >= ref_len || query_index >= seq_len) break;
                
                if (ref[ref_index] == seq[query_index]) {
                    match_count++;
                    if (seq[query_index] == b1) {
                        zf++;
                    } else if (seq[query_index] == b3) {
                        zc++;
                    }
                } else {
                    // Mismatch: add match_count and mismatch base
                    md_pos += sprintf(md_buffer + md_pos, "%d%c", match_count, ref[ref_index]);
                    match_count = 0;
                    
                    if (seq[query_index] == b2) {
                        yf++;
                    } else if (seq[query_index] == b4) {
                        yc++;
                    } else {
                        ns++;
                    }
                }
                ref_index++;
                query_index++;
            }
        } else if (operation == 1) {  // Insertion
            query_index += length;
            nc += length;
        } else if (operation == 4) {  // Soft clipping
            query_index += length;
            nc += length;
        } else if (operation == 2) {  // Deletion
            md_pos += sprintf(md_buffer + md_pos, "%d^", match_count);
            for (long j = 0; j < length && ref_index < ref_len; j++) {
                md_buffer[md_pos++] = ref[ref_index++];
            }
            match_count = 0;
            nc += length;
        }
    }
    
    // Append final match count
    md_pos += sprintf(md_buffer + md_pos, "%d", match_count);
    
    // Create Python string from buffer
    PyObject* md_tag = PyUnicode_FromStringAndSize(md_buffer, md_pos);
    free(md_buffer);
    
    if (!md_tag) {
        return NULL;
    }
    
    // Return tuple (md_tag, yf, zf, yc, zc, ns, nc)
    return Py_BuildValue("(Oiiiiii)", md_tag, yf, zf, yc, zc, ns, nc);
}

// Fast directional score calculation
static PyObject* fast_calculate_directional_score(PyObject* self, PyObject* args) {
    PyObject* cigar_list;
    const char* seq;
    const char* ref;
    int is_orientation1;
    
    if (!PyArg_ParseTuple(args, "Ossp", &cigar_list, &seq, &ref, &is_orientation1)) {
        return NULL;
    }
    
    Py_ssize_t cigar_len = PyList_Size(cigar_list);
    Py_ssize_t seq_len = strlen(seq);
    Py_ssize_t ref_len = strlen(ref);
    
    int ref_index = 0, query_index = 0;
    int matches = 0, expected_conversions = 0, wrong_conversions = 0;
    int other_mismatches = 0, indels = 0;
    
    // Process CIGAR operations
    for (Py_ssize_t i = 0; i < cigar_len; i++) {
        PyObject* cigar_op = PyList_GetItem(cigar_list, i);
        long length, operation;
        
        // Handle both list and tuple formats
        if (PyList_Check(cigar_op) && PyList_Size(cigar_op) == 2) {
            length = PyLong_AsLong(PyList_GetItem(cigar_op, 0));
            operation = PyLong_AsLong(PyList_GetItem(cigar_op, 1));
        } else if (PyTuple_Check(cigar_op) && PyTuple_Size(cigar_op) == 2) {
            length = PyLong_AsLong(PyTuple_GetItem(cigar_op, 0));
            operation = PyLong_AsLong(PyTuple_GetItem(cigar_op, 1));
        } else {
            PyErr_SetString(PyExc_ValueError, "CIGAR must be list of [length, operation] or (length, operation)");
            return NULL;
        }
        
        if (operation == 0) {  // Match or Mismatch
            for (long j = 0; j < length; j++) {
                if (ref_index >= ref_len || query_index >= seq_len) break;
                
                char ref_base = ref[ref_index];
                char read_base = seq[query_index];
                
                if (ref_base == read_base) {
                    matches++;
                } else {
                    // Classify mismatch
                    if (is_orientation1) {
                        // Orientation 1: MK conversion (C->T, A->G expected)
                        if ((ref_base == 'C' && read_base == 'T') || 
                            (ref_base == 'A' && read_base == 'G')) {
                            expected_conversions++;
                        } else if ((ref_base == 'T' && read_base == 'C') || 
                                   (ref_base == 'G' && read_base == 'A')) {
                            wrong_conversions++;
                        } else {
                            other_mismatches++;
                        }
                    } else {
                        // Orientation 2: KM conversion (G->A, T->C expected)
                        if ((ref_base == 'G' && read_base == 'A') || 
                            (ref_base == 'T' && read_base == 'C')) {
                            expected_conversions++;
                        } else if ((ref_base == 'A' && read_base == 'G') || 
                                   (ref_base == 'C' && read_base == 'T')) {
                            wrong_conversions++;
                        } else {
                            other_mismatches++;
                        }
                    }
                }
                ref_index++;
                query_index++;
            }
        } else if (operation == 1 || operation == 2) {  // Insertion or Deletion
            if (operation == 1) {
                query_index += length;
            } else {
                ref_index += length;
            }
            indels += length;
        } else if (operation == 4) {  // Soft clipping
            query_index += length;
        }
    }
    
    // Calculate score
    int score = matches + expected_conversions - wrong_conversions - other_mismatches - indels;
    int total_bad_mismatches = wrong_conversions + other_mismatches;
    
    // Return tuple (score, wrong_conversions, total_bad_mismatches)
    return Py_BuildValue("(iii)", score, wrong_conversions, total_bad_mismatches);
}

// Method definitions
static PyMethodDef SeqOpsMethods[] = {
    {"fast_base_conversion", fast_base_conversion, METH_VARARGS, "Fast base conversion"},
    {"fast_reverse_complement", fast_reverse_complement, METH_VARARGS, "Fast reverse complement"},
    {"fast_cal_md_and_tag", fast_cal_md_and_tag, METH_VARARGS, "Fast MD tag and conversion stats calculation"},
    {"fast_calculate_directional_score", fast_calculate_directional_score, METH_VARARGS, "Fast directional score calculation"},
    {NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef seqopsmodule = {
    PyModuleDef_HEAD_INIT,
    "seqops",
    "Fast sequence operations in C",
    -1,
    SeqOpsMethods
};

// Module initialization
PyMODINIT_FUNC PyInit_seqops(void) {
    return PyModule_Create(&seqopsmodule);
}
