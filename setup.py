from distutils.core import setup

setup(
    name='transitfinder',
    version='0.1.0',
    author='Paul Hildebrandt',
    author_email='me@paulhldbrndt.com',
    url='http://github.com/paulhldbrndt/transitfinder',
    description='transitfinder',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy",
        "astropy",
        "matplotlib",
        "tkinter",
    ],
)
