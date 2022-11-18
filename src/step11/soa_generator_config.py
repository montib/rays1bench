import sys

sys.path.append('../../soa_generator')

import soa_generator as soa


header_extra = '#include "mymath.h"\n\nclass Material;\n'
add_parameters = 'Vec3 center, float radius, Material *material'


properties = [['float', 'center_x', ' = center.getX();'],
                ['float', 'center_y', ' = center.getY();'],
                ['float', 'center_z', ' = center.getZ();'],
                ['float', 'radius_sq', ' = radius*radius;'],
                ['float', 'inv_radius', ' = radius > 0 ? (1.0f / radius) : 0;'],
                ['float', 'discriminant', ' = 0;'],
                ['float', 'nb', ' = 0;'],
                ['Material *', 'material', ' = material;']]

soa.generate('SphereSOA', properties, add_parameters, header_extra, 'soa_sphere')