#!/usr/bin/env python3

from distutils.core import setup
from ancseq.__init__ import __version__

setup(name='ancseq',
      version='{}'.format(__version__),
      description='ancseq: 조상 서열 재구축을 위한 IQ-TREE 래퍼',
      author='Yu Sugihara',
      author_email='sugihara.yu.85s@kyoto-u.jp',
      url='https://github.com/YuSugihara/ancseq',
      license='GPL',
      packages=['ancseq'],
      entry_points={'console_scripts': [
            'ancseq = ancseq.ancseq:main',
            ]
        }
    )
