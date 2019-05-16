import subprocess
from codecs import open
from setuptools import setup, find_packages
from setuptools.command import develop, build_py


def readme():
    with open("README.md", "r", "utf-8") as f:
        return f.read()


class CustomDevelop(develop.develop, object):
    """
    Class needed for "pip install -e ."
    """
    def run(self):
        subprocess.check_call("make", shell=True)
        super(CustomDevelop, self).run()


class CustomBuildPy(build_py.build_py, object):
    """
    Class needed for "pip install s2p"
    """
    def run(self):
        super(CustomBuildPy, self).run()
        subprocess.check_call("make", shell=True)
        subprocess.check_call("cp -r bin lib build/lib/", shell=True)


requirements = ['numpy',
                'scipy',
                'rasterio[s3,test]',
                'utm',
                'pyproj',
                'beautifulsoup4[lxml]',
                'requests']

setup(name="s2p",
      version="1.0b10",
      description="Satellite Stereo Pipeline.",
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://github.com/cmla/s2p',
      packages=['s2p'],
      install_requires=requirements,
      cmdclass={'develop': CustomDevelop,
                'build_py': CustomBuildPy},
      entry_points="""
          [console_scripts]
          s2p=s2p.cli:main
      """)
