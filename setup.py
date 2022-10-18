from setuptools import setup

setup(name='proteinproductionsim',
      version='0.0.1',
      description='TBD',
      url='https://github.com/Zack-Code-UIUC/proteinproductionsim',
      author='Zelong Xiong',
      author_email='zelongx2@illinois.edu',
      license='MIT',
      packages=['proteinproductionsim',
                'proteinproductionsim.controller',
                'proteinproductionsim.datacontainer',
                'proteinproductionsim.entity',
                'proteinproductionsim.environment',
                'proteinproductionsim.helper'
                ],
      zip_safe=False
      )
