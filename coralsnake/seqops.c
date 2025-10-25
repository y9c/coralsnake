#include <Python.h>
#include <string.h>
#include <zlib.h>

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

// Fast FASTA file conversion (line-by-line, memory efficient)
static PyObject* fast_convert_fasta_file(PyObject* self, PyObject* args) {
    const char* input_path;
    const char* output_path;
    const char* from_bases;
    const char* to_bases;
    
    if (!PyArg_ParseTuple(args, "ssss", &input_path, &output_path, &from_bases, &to_bases)) {
        return NULL;
    }
    
    // Open input and output files
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
    
    // Buffer for reading lines (large enough for typical FASTA lines)
    char* line_buffer = (char*)malloc(1048576);  // 1MB buffer
    if (!line_buffer) {
        fclose(input_file);
        fclose(output_file);
        return PyErr_NoMemory();
    }
    
    // Process file line by line
    while (fgets(line_buffer, 1048576, input_file)) {
        if (line_buffer[0] == '>') {
            // Header line - write directly
            fputs(line_buffer, output_file);
        } else {
            // Sequence line - convert bases
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
    
    // Clean up
    free(line_buffer);
    fclose(input_file);
    fclose(output_file);
    
    Py_RETURN_NONE;
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

// Fast FASTQ file splitting for paired-end (supports .fq and .fq.gz)
// Splits R1 and R2 together at exactly the same read positions
static PyObject* fast_split_fastq_paired(PyObject* self, PyObject* args) {
    const char* r1_path;
    const char* r2_path;  // Can be NULL for single-end
    const char* output_dir;
    int num_chunks;
    
    if (!PyArg_ParseTuple(args, "szsi", &r1_path, &r2_path, &output_dir, &num_chunks)) {
        return NULL;
    }
    
    // Check if R1 is gzipped
    int r1_is_gzipped = 0;
    size_t r1_path_len = strlen(r1_path);
    if (r1_path_len > 3 && strcmp(r1_path + r1_path_len - 3, ".gz") == 0) {
        r1_is_gzipped = 1;
    }
    
    // Check if R2 is gzipped (if provided)
    int r2_is_gzipped = 0;
    if (r2_path) {
        size_t r2_path_len = strlen(r2_path);
        if (r2_path_len > 3 && strcmp(r2_path + r2_path_len - 3, ".gz") == 0) {
            r2_is_gzipped = 1;
        }
    }
    
    // Get R1 file size to estimate number of reads
    FILE* size_file = fopen(r1_path, "rb");
    if (!size_file) {
        PyErr_SetString(PyExc_IOError, "Cannot open R1 file");
        return NULL;
    }
    fseek(size_file, 0, SEEK_END);
    long file_size = ftell(size_file);
    fclose(size_file);
    
    // Estimate reads per chunk based on R1 file size
    // Assume average read length ~150bp, with quality and headers ~400 bytes per read
    // For gzipped files, assume 3-4x compression
    long estimated_reads = r1_is_gzipped ? (file_size * 3 / 400) : (file_size / 400);
    long reads_per_chunk = (estimated_reads + num_chunks - 1) / num_chunks;
    
    // Buffers for reading lines
    char* r1_line_buffer = (char*)malloc(65536);  // 64KB buffer for R1
    char* r2_line_buffer = r2_path ? (char*)malloc(65536) : NULL;  // 64KB buffer for R2 if needed
    if (!r1_line_buffer || (r2_path && !r2_line_buffer)) {
        free(r1_line_buffer);
        free(r2_line_buffer);
        return PyErr_NoMemory();
    }
    
    // Create output file handles for R1
    FILE** r1_output_files = (FILE**)malloc(num_chunks * sizeof(FILE*));
    FILE** r2_output_files = r2_path ? (FILE**)malloc(num_chunks * sizeof(FILE*)) : NULL;
    
    if (!r1_output_files || (r2_path && !r2_output_files)) {
        free(r1_line_buffer);
        free(r2_line_buffer);
        free(r1_output_files);
        free(r2_output_files);
        return PyErr_NoMemory();
    }
    
    // Open all R1 output files
    for (int i = 0; i < num_chunks; i++) {
        char output_path[1024];
        snprintf(output_path, sizeof(output_path), "%s/chunk_%d_R1.fq", output_dir, i);
        r1_output_files[i] = fopen(output_path, "w");
        if (!r1_output_files[i]) {
            for (int j = 0; j < i; j++) {
                fclose(r1_output_files[j]);
            }
            free(r1_output_files);
            free(r2_output_files);
            free(r1_line_buffer);
            free(r2_line_buffer);
            PyErr_SetString(PyExc_IOError, "Cannot open R1 output file");
            return NULL;
        }
    }
    
    // Open all R2 output files if needed
    if (r2_path) {
        for (int i = 0; i < num_chunks; i++) {
            char output_path[1024];
            snprintf(output_path, sizeof(output_path), "%s/chunk_%d_R2.fq", output_dir, i);
            r2_output_files[i] = fopen(output_path, "w");
            if (!r2_output_files[i]) {
                for (int j = 0; j < i; j++) {
                    fclose(r2_output_files[j]);
                }
                for (int j = 0; j < num_chunks; j++) {
                    fclose(r1_output_files[j]);
                }
                free(r1_output_files);
                free(r2_output_files);
                free(r1_line_buffer);
                free(r2_line_buffer);
                PyErr_SetString(PyExc_IOError, "Cannot open R2 output file");
                return NULL;
            }
        }
    }
    
    // Open R1 file
    void* r1_file;
    if (r1_is_gzipped) {
        r1_file = gzopen(r1_path, "rb");
        if (!r1_file) {
            // Cleanup and error
            for (int i = 0; i < num_chunks; i++) {
                fclose(r1_output_files[i]);
                if (r2_path) fclose(r2_output_files[i]);
            }
            free(r1_output_files);
            free(r2_output_files);
            free(r1_line_buffer);
            free(r2_line_buffer);
            PyErr_SetString(PyExc_IOError, "Cannot open R1 file");
            return NULL;
        }
    } else {
        r1_file = fopen(r1_path, "r");
        if (!r1_file) {
            // Cleanup and error
            for (int i = 0; i < num_chunks; i++) {
                fclose(r1_output_files[i]);
                if (r2_path) fclose(r2_output_files[i]);
            }
            free(r1_output_files);
            free(r2_output_files);
            free(r1_line_buffer);
            free(r2_line_buffer);
            PyErr_SetString(PyExc_IOError, "Cannot open R1 file");
            return NULL;
        }
    }
    
    // Open R2 file if provided
    void* r2_file = NULL;
    if (r2_path) {
        if (r2_is_gzipped) {
            r2_file = gzopen(r2_path, "rb");
        } else {
            r2_file = fopen(r2_path, "r");
        }
        if (!r2_file) {
            // Cleanup and error
            if (r1_is_gzipped) gzclose(r1_file);
            else fclose(r1_file);
            for (int i = 0; i < num_chunks; i++) {
                fclose(r1_output_files[i]);
                fclose(r2_output_files[i]);
            }
            free(r1_output_files);
            free(r2_output_files);
            free(r1_line_buffer);
            free(r2_line_buffer);
            PyErr_SetString(PyExc_IOError, "Cannot open R2 file");
            return NULL;
        }
    }
    
    // Read and distribute reads (R1 and R2 together)
    int current_chunk = 0;
    long reads_in_current_chunk = 0;
    int line_in_read = 0;
    
    while (1) {
        // Read one line from R1
        char* r1_line;
        if (r1_is_gzipped) {
            r1_line = gzgets((gzFile)r1_file, r1_line_buffer, 65536);
        } else {
            r1_line = fgets(r1_line_buffer, 65536, (FILE*)r1_file);
        }
        
        if (!r1_line) break;  // End of R1 file
        
        // Write R1 line to current chunk
        fputs(r1_line_buffer, r1_output_files[current_chunk]);
        
        // Read and write corresponding R2 line if paired-end
        if (r2_path) {
            char* r2_line;
            if (r2_is_gzipped) {
                r2_line = gzgets((gzFile)r2_file, r2_line_buffer, 65536);
            } else {
                r2_line = fgets(r2_line_buffer, 65536, (FILE*)r2_file);
            }
            
            if (r2_line) {
                fputs(r2_line_buffer, r2_output_files[current_chunk]);
            }
        }
        
        line_in_read++;
        if (line_in_read == 4) {
            // Completed one read (4 lines)
            line_in_read = 0;
            reads_in_current_chunk++;
            
            // Move to next chunk if current is full (except for last chunk)
            if (reads_in_current_chunk >= reads_per_chunk && current_chunk < num_chunks - 1) {
                current_chunk++;
                reads_in_current_chunk = 0;
            }
        }
    }
    
    // Close input files
    if (r1_is_gzipped) gzclose((gzFile)r1_file);
    else fclose((FILE*)r1_file);
    
    if (r2_path) {
        if (r2_is_gzipped) gzclose((gzFile)r2_file);
        else fclose((FILE*)r2_file);
    }
    
    // Close all output files
    for (int i = 0; i < num_chunks; i++) {
        fclose(r1_output_files[i]);
        if (r2_path) fclose(r2_output_files[i]);
    }
    
    // Cleanup
    free(r1_output_files);
    free(r2_output_files);
    free(r1_line_buffer);
    free(r2_line_buffer);
    
    // Return list of tuples: [(r1_chunk0, r2_chunk0), ...] or [(r1_chunk0, None), ...] for single-end
    PyObject* result = PyList_New(num_chunks);
    for (int i = 0; i < num_chunks; i++) {
        char r1_path_out[1024];
        snprintf(r1_path_out, sizeof(r1_path_out), "%s/chunk_%d_R1.fq", output_dir, i);
        
        if (r2_path) {
            char r2_path_out[1024];
            snprintf(r2_path_out, sizeof(r2_path_out), "%s/chunk_%d_R2.fq", output_dir, i);
            PyObject* tuple = PyTuple_Pack(2, 
                PyUnicode_FromString(r1_path_out),
                PyUnicode_FromString(r2_path_out));
            PyList_SetItem(result, i, tuple);
        } else {
            PyObject* tuple = PyTuple_Pack(2,
                PyUnicode_FromString(r1_path_out),
                Py_None);
            Py_INCREF(Py_None);
            PyList_SetItem(result, i, tuple);
        }
    }
    
    return result;
}

// Method definitions
static PyMethodDef SeqOpsMethods[] = {
    {"fast_base_conversion", fast_base_conversion, METH_VARARGS, "Fast base conversion"},
    {"fast_reverse_complement", fast_reverse_complement, METH_VARARGS, "Fast reverse complement"},
    {"fast_cal_md_and_tag", fast_cal_md_and_tag, METH_VARARGS, "Fast MD tag and conversion stats calculation"},
    {"fast_calculate_directional_score", fast_calculate_directional_score, METH_VARARGS, "Fast directional score calculation"},
    {"fast_convert_fasta_file", fast_convert_fasta_file, METH_VARARGS, "Fast FASTA file conversion (line-by-line)"},
    {"fast_split_fastq_paired", fast_split_fastq_paired, METH_VARARGS, "Fast paired-end FASTQ file splitting (supports .fq and .fq.gz)"},
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
