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

// Fast FASTQ file splitting (supports .fq and .fq.gz)
static PyObject* fast_split_fastq(PyObject* self, PyObject* args) {
    const char* input_path;
    const char* output_dir;
    int num_chunks;
    int is_r2;  // 0 for R1, 1 for R2
    
    if (!PyArg_ParseTuple(args, "ssip", &input_path, &output_dir, &num_chunks, &is_r2)) {
        return NULL;
    }
    
    // Check if gzipped by file extension
    int is_gzipped = 0;
    size_t path_len = strlen(input_path);
    if (path_len > 3 && strcmp(input_path + path_len - 3, ".gz") == 0) {
        is_gzipped = 1;
    }
    
    // Get file size to estimate number of reads
    FILE* size_file = fopen(input_path, "rb");
    if (!size_file) {
        PyErr_SetString(PyExc_IOError, "Cannot open input file");
        return NULL;
    }
    fseek(size_file, 0, SEEK_END);
    long file_size = ftell(size_file);
    fclose(size_file);
    
    // Estimate reads per chunk
    // Assume average read length ~150bp, with quality and headers ~400 bytes per read
    // For gzipped files, assume 3-4x compression
    long estimated_reads = is_gzipped ? (file_size * 3 / 400) : (file_size / 400);
    long reads_per_chunk = (estimated_reads + num_chunks - 1) / num_chunks;
    
    // Buffer for reading lines
    char* line_buffer = (char*)malloc(65536);  // 64KB buffer
    if (!line_buffer) {
        return PyErr_NoMemory();
    }
    
    // Create output file handles
    FILE** output_files = (FILE**)malloc(num_chunks * sizeof(FILE*));
    if (!output_files) {
        free(line_buffer);
        return PyErr_NoMemory();
    }
    
    // Open all output files
    const char* suffix = is_r2 ? "_R2.fq" : "_R1.fq";
    for (int i = 0; i < num_chunks; i++) {
        char output_path[1024];
        snprintf(output_path, sizeof(output_path), "%s/chunk_%d%s", output_dir, i, suffix);
        output_files[i] = fopen(output_path, "w");
        if (!output_files[i]) {
            // Clean up opened files
            for (int j = 0; j < i; j++) {
                fclose(output_files[j]);
            }
            free(output_files);
            free(line_buffer);
            PyErr_SetString(PyExc_IOError, "Cannot open output file");
            return NULL;
        }
    }
    
    // Read and distribute reads
    int current_chunk = 0;
    long reads_in_current_chunk = 0;
    int line_in_read = 0;  // 0-3 for the 4 lines of a FASTQ read
    
    if (is_gzipped) {
        // Use zlib for gzipped files
        gzFile gz_file = gzopen(input_path, "rb");
        if (!gz_file) {
            for (int i = 0; i < num_chunks; i++) {
                fclose(output_files[i]);
            }
            free(output_files);
            free(line_buffer);
            PyErr_SetString(PyExc_IOError, "Cannot open gzipped input file");
            return NULL;
        }
        
        while (gzgets(gz_file, line_buffer, 65536)) {
            // Write line to current chunk
            fputs(line_buffer, output_files[current_chunk]);
            
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
        gzclose(gz_file);
    } else {
        // Use regular fopen for uncompressed files
        FILE* input_file = fopen(input_path, "r");
        if (!input_file) {
            for (int i = 0; i < num_chunks; i++) {
                fclose(output_files[i]);
            }
            free(output_files);
            free(line_buffer);
            PyErr_SetString(PyExc_IOError, "Cannot open input file");
            return NULL;
        }
        
        while (fgets(line_buffer, 65536, input_file)) {
            // Write line to current chunk
            fputs(line_buffer, output_files[current_chunk]);
            
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
        fclose(input_file);
    }
    
    // Close all files
    for (int i = 0; i < num_chunks; i++) {
        fclose(output_files[i]);
    }
    free(output_files);
    free(line_buffer);
    
    // Return list of output file paths
    PyObject* result = PyList_New(num_chunks);
    for (int i = 0; i < num_chunks; i++) {
        char output_path[1024];
        snprintf(output_path, sizeof(output_path), "%s/chunk_%d%s", output_dir, i, suffix);
        PyList_SetItem(result, i, PyUnicode_FromString(output_path));
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
    {"fast_split_fastq", fast_split_fastq, METH_VARARGS, "Fast FASTQ file splitting (supports .fq and .fq.gz)"},
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
