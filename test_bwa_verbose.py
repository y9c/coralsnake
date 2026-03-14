from bwamem.libbwa import libbwa, ffi

try:
    libbwa.bwa_verbose = 1
    print("Set bwa_verbose successfully!")
except Exception as e:
    print(e)
