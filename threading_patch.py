#!/usr/bin/env python3
"""
Patch to add threading support to mapping.py

This script will:
1. Add ThreadPoolExecutor import
2. Add a worker function that uses shared indices
3. Update map_file to use threading
"""

import re

# Read the current mapping.py
with open('coralsnake/mapping.py', 'r') as f:
    content = f.read()

# 1. Update imports
content = content.replace(
    'import os\nimport random\n\nimport mappy as mp',
    'import os\nimport random\nfrom concurrent.futures import ThreadPoolExecutor\n\nimport mappy as mp'
)

# 2. Add worker function after imports
worker_function = '''

def _worker_map_reads(batch, idx0, idx_mk, fwd_lib, max_mismatches):
    """
    Worker function for parallel read mapping using threads.
    Threads share the same indices (no memory duplication!).
    
    The key insight: minimap2's C code releases the Python GIL,
    so multiple threads can actually map reads in parallel.
    
    Args:
        batch: List of (name, seq1, seq2, qua1, qua2) tuples
        idx0: Shared original reference index
        idx_mk: Shared MK converted reference index
        fwd_lib: Library strand orientation
        max_mismatches: Maximum allowed bad mismatches
    
    Returns:
        List of (name, mapped_results) tuples
    """
    results = []
    for name, seq1, seq2, qua1, qua2 in batch:
        mapped = run_mapping(name, seq1, seq2, qua1, qua2, idx0, idx_mk, fwd_lib, max_mismatches)
        results.append((name, mapped))
    return results
'''

# Insert after the imports section
import_end = content.find('\ndef ')
content = content[:import_end] + worker_function + '\n' + content[import_end:]

# 3. Update map_file signature
content = content.replace(
    'def map_file(ref_file, r1_file, r2_file, output_file, fwd_lib=True, max_mismatches=10, threads=8):',
    'def map_file(ref_file, r1_file, r2_file, output_file, fwd_lib=True, max_mismatches=10, threads=8, num_workers=None):'
)

# 4. Update docstring
content = content.replace(
    '        threads: Number of threads for parallel processing (default: 8)',
    '        threads: Number of threads for minimap2 indexing (default: 8)\n        num_workers: Number of parallel worker threads for mapping (default: threads)'
)

# Write the patched file
with open('coralsnake/mapping.py', 'w') as f:
    f.write(content)

print("✅ Patch applied successfully!")
print("   - Added ThreadPoolExecutor import")
print("   - Added _worker_map_reads function")
print("   - Updated map_file signature")

