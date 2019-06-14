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


try:
    from wheel.bdist_wheel import bdist_wheel
    class BdistWheel(bdist_wheel):
        """
        Class needed to build platform dependent binary wheels
        """
        def finalize_options(self):
            bdist_wheel.finalize_options(self)
            self.root_is_pure = False

except ImportError:
    BdistWheel = None


requirements = ['numpy',
                'scipy',
                'urllib3<1.25',  # This is a temporary hack to fix https://github.com/boto/botocore/issues/1733 until next release of botocore
                'rasterio[s3,test]>=1.0.14',
                'utm',
                'pyproj',
                'beautifulsoup4[lxml]',
                'plyfile',
                'ransac',
                'requests']

setup(name="s2p",
      version="1.0b16",
      description="Satellite Stereo Pipeline.",
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://github.com/cmla/s2p',
      packages=['s2p'],
      install_requires=requirements,
      cmdclass={'develop': CustomDevelop,
                'build_py': CustomBuildPy,
                'bdist_wheel': BdistWheel},
      entry_points="""
          [console_scripts]
          s2p=s2p.cli:main
      """)
