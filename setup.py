from setuptools import setup


def readme():
    with open('README.md', 'r') as f:
        content = f.read()
    return content

setup(
    name = "bmk-ocean-stars",
    version = 'v0.6.9-beta',
    packages = ["bmkos"],
    python_requires='>=3.8, <3.10',
    install_requires = ['numpy', 'scipy', 'pandas', 'tqdm', 'setuptools >= 18.0',
                        'editdistance', 'pysam', 'parasail', 'bioframe >= 0.3.0',
                        'jinja2', 'plotly', 'bmk-cell-calling >= 0.2.1'],

    author = "Wang Jianshou",
    author_email = "wangjs@biomarker.com.cn",
    description='Analysis of nanopore sequencing data of single cell.',
    long_description=readme(),
    long_description_content_type='text/markdown',
    include_package_data=True,
    url='https://github.com/wangjianshou/bmk-ocean-stars',
    entry_points={
              'console_scripts': [
                  'bmkos = bmkos.__main__:main',
              ]
        },
    classifiers=[
              'Development Status :: 3 - Alpha',
              'Operating System :: POSIX :: Linux',
              'Programming Language :: Python :: 3.9',
              'Programming Language :: Python :: 3.8',
              #'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        ],
)
