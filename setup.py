import subprocess
from setuptools import setup, find_packages
from setuptools.command import develop, build_py


def readme():
    with open('README.md') as f:
        return f.read()


class CustomDevelop(develop.develop):
    """
    Class needed for "pip install -e ."
    """
    def run(self):
        subprocess.check_call("make CC=gcc-8 CXX=g++-8", shell=True)
        super().run()


class CustomBuildPy(build_py.build_py):
    """
    Class needed for "pip install pys2p"
    """
    def run(self):
        super().run()
        subprocess.check_call("make", shell=True)
        subprocess.check_call("cp -r bin build/lib/", shell=True)


requirements = ['numpy',
                'scipy',
                'rasterio[s3]',
                'utm',
                'pyproj',
                'bs4',
                'requests']

setup(name="s2p",
      version="1.0a1",
      description="Satellite Stereo Pipeline.",
      long_description=readme(),
      url='https://github.com/miss3d/s2p',
      packages=['s2p'],
      install_requires=requirements,
      cmdclass={'develop': CustomDevelop,
                'build_py': CustomBuildPy},
      entry_points="""
          [console_scripts]
          s2p=s2p.cli:main
      """)
