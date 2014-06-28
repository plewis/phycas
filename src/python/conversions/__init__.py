from _ConversionsExt import *

import os
from phycas import release_version
print 'release_version is',release_version
if not release_version:
    print 'importing conversions from ',os.path.abspath(__file__)
