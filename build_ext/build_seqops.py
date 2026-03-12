from setuptools import setup, Extension

setup(
    name="seqops",
    ext_modules=[
        Extension(
            "coralsnake.seqops",
            ["coralsnake/seqops.c"],
            libraries=["z"],
        )
    ],
)
