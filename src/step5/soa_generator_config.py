import sys

sys.path.append('../../soa_generator')

import soa_generator as soa


header_extra = '#include "mymath.h"\n\nclass Material;\n'
add_parameters = 'Vec3 center, float radius, Material *material'

properties = [['Vec3', 'center', ' = center;'],
                ['float', 'radius', ' = radius;'],
                ['Material *', 'material', ' = material;']]

soa.generate('SphereSOA', properties, add_parameters, header_extra, 'soa_sphere')