from string import Template
import os

kType = 0
kName = 1
kDefualtValue = 2

indent = '    '


def generate_header(class_name, properties, add_parameters_str, header_extra):
    template = """\
#pragma once

// generated file!

#include <stdint.h>
#include <assert.h>

${HeaderExtra}


class ${ClassName}
{
public:
    ${ClassName}();
    ~${ClassName}();

    void cleanUp();

    uint32_t add(${AddParameterList});
    void copy(uint32_t dst_idx, uint32_t src_id);
    void remove(uint32_t idx);                      // ! swap with last and pop
    void reserve(uint32_t capacity);

    inline uint32_t getCount() const
    {
        return _data._count;
    }

    inline void setCount(uint32_t new_size)         // decrement allowed only; use add() to add an element
    {
        assert(new_size <= _data._count);
        _data._count = new_size;
    }

    struct InstanceData
    {
        ${InstanceDataSize}

        void *_buffer;
        uint32_t _count;
        uint32_t _capacity;

        ${InstanceDataMembers}
    };

    InstanceData *getData()
    {
        return &_data;
    }

    const InstanceData *getData() const
    {
        return &_data;
    }

private:
    void allocate(uint32_t new_count);

    InstanceData _data;
};
"""

    prop_str = ''
    instance_data_size_str = 'enum { Size = '

    for prop in properties:
        name = prop[kName]
        type_name = prop[kType]

        prop_str += '\n' + indent * 2 + type_name + ('*' if type_name[-1] == '*' else ' *') + ' __restrict ' + name + ';'

        instance_data_size_str += ' + sizeof(' + type_name + ')'

    instance_data_size_str += ' };'

    return Template(template).substitute(
        ClassName=class_name,
        AddParameterList=add_parameters_str,
        InstanceDataSize=instance_data_size_str,
        InstanceDataMembers=prop_str,
        HeaderExtra=header_extra)


def generate_cpp(class_name, properties, add_parameters_str, filename_h):
    template = """\
// generated file!

#include "${HeaderFileName}"
#include "allocator.h"
#include <string.h>


static inline bool isPower2(size_t x)
{
    return x && !(x & (x - 1));
}

static inline char *align_ptr(char *ptr, size_t align)
{
    assert(align > 0 && isPower2((uint32_t)align));

    return (char*)(((uintptr_t)ptr + (align-1)) & ~(uintptr_t)(align-1));
}



${ClassName}::${ClassName}()
{
    memset(&_data, 0, sizeof(_data));
}


${ClassName}::~${ClassName}()
{
    cleanUp();
}

void ${ClassName}::cleanUp()
{
    if (_data._buffer)
    {
        aligned_free(_data._buffer);
        _data._buffer = nullptr;${CleanUpPointers}
    }

    memset(&_data, 0, sizeof(_data));
}

void ${ClassName}::copy(uint32_t dst_idx, uint32_t src_idx)
{${CopyBlock}
}

void ${ClassName}::remove(uint32_t idx)
{
    --_data._count;
    uint32_t last_idx = _data._count;
    uint32_t del_idx = idx;

    if (del_idx != last_idx)
        copy(del_idx, last_idx);
}

uint32_t ${ClassName}::add(${AddParameterList})
{
    reserve(_data._count+1);

    uint32_t idx = _data._count++;

    // init
${AddDefaults}
    return idx;
}

void ${ClassName}::reserve(uint32_t capacity)
{
    if (_data._capacity < capacity)
        allocate(1 + 2*capacity);
}

void ${ClassName}::allocate(uint32_t new_capacity)
{
    assert(new_capacity > _data._count);    // it's a simple container :)

    const size_t Align = 32;
    InstanceData new_data;
    ${AllocateDefineNewBytes}

    new_data._buffer = aligned_allocate(new_bytes, Align);
    new_data._capacity = new_capacity;
    new_data._count = _data._count;
    ${AllocateSetPointers}

    if (_data._count)
    {${AllocateCopy}
    }

    if (_data._buffer)
       aligned_free(_data._buffer);

    _data = new_data;
}

"""
    prev_property_name = ''
    prev_type_name = ''
    cleanup_pointers_str = ''
    allocate_set_pointers_str = ''
    allocate_copy_str = ''
    add_defaults_str = ''

    copy_block_str = ''

    for i in range(len(properties)):
        prop = properties[i]
        name = prop[kName]
        type_name = prop[kType]

        cleanup_pointers_str += '\n' + indent * 2 + '_data.' + name + ' = nullptr;'
        copy_block_str += '\n' + indent + '_data.' + name + '[dst_idx] = _data.' + name + '[src_idx];'

        if i > 0:
            prev_property_name = properties[i - 1][kName]
            prev_type_name = properties[i - 1][kType]
            allocate_set_pointers_str += '\n' + indent + 'new_data.' + name + ' = (' + type_name + '*)(align_ptr((char*)&new_data.' + prev_property_name + '[new_capacity], Align));'
        else:
            allocate_set_pointers_str += '\n' + indent + 'new_data.' + name + ' = (' + type_name + '*)new_data._buffer;'

        # # asserts
        # allocate_copy_str += '\n\n' + indent * 2 + 'assert((uintptr_t)new_data.' + name + ' % alignof(' + type_name + ') == 0);'

        # if i > 0:
        #     allocate_copy_str += '\n' + indent * 2 + 'assert((uintptr_t)new_data.' + name + ' - (uintptr_t)&new_data.' + prev_property_name + '[new_capacity] <= Align);'
        #     allocate_copy_str += '\n' + indent * 2 + 'assert((uintptr_t)new_data.' + name + ' >= (uintptr_t)&new_data.' + prev_property_name + '[new_capacity]);'
        #     allocate_copy_str += '\n' + indent * 2 + 'assert((uintptr_t)new_data.' + name + ' > new_capacity * sizeof(' + prev_type_name + '));'

        # copy data
        allocate_copy_str += '\n' + indent * 2 + 'memcpy(new_data.' + name + ', _data.' + name + ', _data._count * sizeof(' + type_name + '));'

        if prop[kDefualtValue] is not None:
            add_defaults_str += indent +  '_data.' + name + '[idx]' + prop[kDefualtValue] + '\n'

    define_new_bytes_str = 'const size_t new_bytes = new_capacity * InstanceData::Size + Align * ' + str(len(properties)) + ';'

    return Template(template).substitute(
        ClassName=class_name,
        HeaderFileName=filename_h,
        CleanUpPointers=cleanup_pointers_str,
        CopyBlock=copy_block_str,
        AllocateDefineNewBytes=define_new_bytes_str,
        AllocateSetPointers=allocate_set_pointers_str,
        AllocateCopy=allocate_copy_str,
        AddDefaults=add_defaults_str,
        AddParameterList=add_parameters_str)


def generate(class_name, properties, add_parameters_str, header_extra, path):
    filename_cpp = path + '.cpp'
    filename_h = path + '.h'

    header = generate_header(class_name, properties, add_parameters_str, header_extra)
    cpp = generate_cpp(class_name, properties, add_parameters_str, os.path.basename(filename_h))

    with open(filename_cpp, 'wt', newline='\n') as fcpp, open(filename_h, 'wt', newline='\n') as fh:
        fcpp.write(cpp)
        fh.write(header)

    print(filename_cpp)
    print(filename_h)
