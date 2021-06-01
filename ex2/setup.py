from setuptools import setup, Extension, find_packages

setup(name='mykmeanssp',
      ext_modules=
        [Extension('mykmeanssp', 
                    ['kmeans.c'],
                    )
        ]
      )
