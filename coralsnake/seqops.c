#include <Python.h>
#include <string.h>
#include <zlib.h>

// Base conversion function
static PyObject* base_conversion(PyObject* self, PyObject* args) {
    const char* seq;
    const char* from_bases;
    const char* to_bases;
    Py_ssize_t seq_len;
    
    if (!PyArg_ParseTuple(args, "sss", &seq, &from_bases, &to_bases)) {
        return NULL;
    }
    
    seq_len = strlen(seq);
    char* target = (char*)malloc(seq_len + 1);
    if (!target) {
        return PyErr_NoMemory();
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
        target[i] = lookup[c] ? lookup[c] : c;
    }
    target[seq_len] = '\0';
    
    PyObject* result = PyUnicode_FromStringAndSize(target, seq_len);
    free(target);
    return result;
}

// Reverse complement function
static PyObject* reverse_complement(PyObject* self, PyObject* args) {
    const char* seq;
    Py_ssize_t seq_len;
    
    if (!PyArg_ParseTuple(args, "s", &seq)) {
        return NULL;
    }
    
    seq_len = strlen(seq);
    char* target = (char*)malloc(seq_len + 1);
    if (!target) {
        return PyErr_NoMemory();
    }
    
    // Complement lookup table
    char complement[256];
    memset(complement, 0, 256);
    complement['A'] = 'T'; complement['T'] = 'A';
    complement['G'] = 'C'; complement['C'] = 'G';
    complement['N'] = 'N';
    complement['a'] = 't'; complement['t'] = 'a';
    complement['g'] = 'c'; complement['c'] = 'g';
    complement['n'] = 'n';
    
    // Reverse and complement
    for (Py_ssize_t i = 0; i < seq_len; i++) {
        unsigned char c = (unsigned char)seq[seq_len - 1 - i];
        target[i] = complement[c] ? complement[c] : c;
    }
    target[seq_len] = '\0';
    
    PyObject* result = PyUnicode_FromStringAndSize(target, seq_len);
    free(target);
    return result;
}

// Optimized combined function: Score + MD tag + Stats
static PyObject* cal_score_and_md_and_tag(PyObject* self, PyObject* args) {
    const char* cigar_str;
    const char* seq;
    const char* ref;
    int is_orientation1;
    
    if (!PyArg_ParseTuple(args, "sssp", &cigar_str, &seq, &ref, &is_orientation1)) {
        return NULL;
    }
    
    Py_ssize_t seq_len = strlen(seq);
    Py_ssize_t ref_len = strlen(ref);
    
    size_t buffer_size = seq_len * 10 + 64; 
    char* md_buffer = (char*)malloc(buffer_size);
    if (!md_buffer) {
        return PyErr_NoMemory();
    }
    
    int yf = 0, zf = 0, yc = 0, zc = 0, ns = 0, nc = 0;
    int matches = 0, expected_conversions = 0, wrong_conversions = 0;
    int other_mismatches = 0, indels = 0;
    int ref_index = 0, query_index = 0, match_count = 0;
    int md_pos = 0;
    
    // Orientation 1: MK conversion (C->T, A->G expected)
    // Orientation 2: KM conversion (G->A, T->C expected)
    char b_match1, b_conv1, b_match2, b_conv2;
    if (is_orientation1) {
        b_match1 = 'A'; b_conv1 = 'G'; // A -> G
        b_match2 = 'C'; b_conv2 = 'T'; // C -> T
    } else {
        b_match1 = 'G'; b_conv1 = 'A'; // G -> A
        b_match2 = 'T'; b_conv2 = 'C'; // T -> C
    }
    
    long length = 0;
    for (int i = 0; cigar_str[i]; i++) {
        if (cigar_str[i] >= '0' && cigar_str[i] <= '9') {
            length = length * 10 + (cigar_str[i] - '0');
        } else {
            char op_char = cigar_str[i];
            int operation = -1;
            switch(op_char) {
                case 'M': operation = 0; break;
                case 'I': operation = 1; break;
                case 'D': operation = 2; break;
                case 'N': operation = 3; break;
                case 'S': operation = 4; break;
                case 'H': operation = 5; break;
                case 'P': operation = 6; break;
                case '=': operation = 7; break;
                case 'X': operation = 8; break;
            }
            
            if (operation == 0 || operation == 7 || operation == 8) {  // M, =, X
                for (long j = 0; j < length; j++) {
                    if (ref_index >= ref_len || query_index >= seq_len) break;
                    
                    char r_b = ref[ref_index];
                    if (r_b >= 'a' && r_b <= 'z') r_b -= 32;
                    char s_b = seq[query_index];
                    if (s_b >= 'a' && s_b <= 'z') s_b -= 32;

                    if (r_b == s_b) {
                        match_count++;
                        matches++;
                        if (s_b == b_match1) {
                            zf++;
                        } else if (s_b == b_match2) {
                            zc++;
                        }
                    } else {
                        // Mismatch: add match_count and mismatch base to MD
                        md_pos += snprintf(md_buffer + md_pos, buffer_size - md_pos, "%d%c", match_count, r_b);
                        match_count = 0;
                        
                        // Classify mismatch for scoring
                        if ((r_b == b_match1 && s_b == b_conv1) || (r_b == b_match2 && s_b == b_conv2)) {
                            expected_conversions++;
                            if (s_b == b_conv1) {
                                yf++;
                            } else {
                                yc++;
                            }
                        } else if ((r_b == b_conv1 && s_b == b_match1) || (r_b == b_conv2 && s_b == b_match2)) {
                            wrong_conversions++;
                            ns++;
                        } else {
                            other_mismatches++;
                            ns++;
                        }
                    }
                    ref_index++;
                    query_index++;
                }
            } else if (operation == 1) {  // Insertion
                query_index += length;
                nc += length;
                indels += length;
            } else if (operation == 4) {  // Soft clipping
                query_index += length;
                nc += length;
            } else if (operation == 2) {  // Deletion
                md_pos += snprintf(md_buffer + md_pos, buffer_size - md_pos, "%d^", match_count);
                for (long j = 0; j < length && ref_index < ref_len; j++) {
                    char r_b = ref[ref_index++];
                    if (r_b >= 'a' && r_b <= 'z') r_b -= 32;
                    if (md_pos < buffer_size - 1) {
                        md_buffer[md_pos++] = r_b;
                    }
                }
                match_count = 0;
                nc += length;
                indels += length;
            } else if (operation == 3) { // N (skip)
                ref_index += length;
            }
            
            length = 0;
        }
    }
    
    // Append final match count to MD
    md_pos += snprintf(md_buffer + md_pos, buffer_size - md_pos, "%d", match_count);
    
    // Calculate final score
    int score = matches + expected_conversions - wrong_conversions - other_mismatches - indels;
    int bad_mm = wrong_conversions + other_mismatches;
    
    PyObject* md_tag = PyUnicode_FromStringAndSize(md_buffer, md_pos);
    free(md_buffer);
    
    if (!md_tag) {
        return NULL;
    }
    
    // Return tuple: (score, bad_mm, md_tag, yf, zf, yc, zc, ns, nc)
    return Py_BuildValue("(iiNiiiiii)", score, bad_mm, md_tag, yf, zf, yc, zc, ns, nc);
}

// FASTA file conversion (line-by-line, memory efficient)
static PyObject* convert_fasta_file(PyObject* self, PyObject* args) {
    const char* input_path;
    const char* output_path;
    const char* from_bases;
    const char* to_bases;
    
    if (!PyArg_ParseTuple(args, "ssss", &input_path, &output_path, &from_bases, &to_bases)) {
        return NULL;
    }
    
    FILE* input_file = fopen(input_path, "r");
    if (!input_file) {
        PyErr_SetString(PyExc_IOError, "Cannot open input file");
        return NULL;
    }
    
    FILE* output_file = fopen(output_path, "w");
    if (!output_file) {
        fclose(input_file);
        PyErr_SetString(PyExc_IOError, "Cannot open output file");
        return NULL;
    }
    
    char lookup[256];
    memset(lookup, 0, 256);
    for (int i = 0; from_bases[i] && to_bases[i]; i++) {
        lookup[(unsigned char)from_bases[i]] = to_bases[i];
        if (from_bases[i] >= 'A' && from_bases[i] <= 'Z') {
            lookup[(unsigned char)from_bases[i] + 32] = to_bases[i] + 32;
        } else if (from_bases[i] >= 'a' && from_bases[i] <= 'z') {
            lookup[(unsigned char)from_bases[i] - 32] = to_bases[i] - 32;
        }
    }
    
    char* line_buffer = (char*)malloc(1048576);  // 1MB buffer
    if (!line_buffer) {
        fclose(input_file);
        fclose(output_file);
        return PyErr_NoMemory();
    }
    
    Py_BEGIN_ALLOW_THREADS
    while (fgets(line_buffer, 1048576, input_file)) {
        if (line_buffer[0] == '>') {
            fputs(line_buffer, output_file);
        } else {
            size_t len = strlen(line_buffer);
            for (size_t i = 0; i < len; i++) {
                unsigned char c = (unsigned char)line_buffer[i];
                if (lookup[c]) {
                    line_buffer[i] = lookup[c];
                }
            }
            fputs(line_buffer, output_file);
        }
    }
    Py_END_ALLOW_THREADS
    
    free(line_buffer);
    fclose(input_file);
    fclose(output_file);
    
    Py_RETURN_NONE;
}

// Method definitions
static PyMethodDef SeqOpsMethods[] = {
    {"base_conversion", base_conversion, METH_VARARGS, "Base conversion"},
    {"reverse_complement", reverse_complement, METH_VARARGS, "Reverse complement"},
    {"cal_score_and_md_and_tag", cal_score_and_md_and_tag, METH_VARARGS, "Combined score, MD tag, and stats calculation"},
    {"convert_fasta_file", convert_fasta_file, METH_VARARGS, "FASTA file conversion (line-by-line)"},
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
