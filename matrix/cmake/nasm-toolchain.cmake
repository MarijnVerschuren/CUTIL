
if (NOT DEFINED CMAKE_ASM_NASM_OBJECT_FORMAT)
  set(CMAKE_ASM_NASM_OBJECT_FORMAT elf64)
endif ()

# nasm output format
if (CMAKE_ASM_NASM_OBJECT_FORMAT MATCHES "elf64")
  set(extra_flags "-m elf_x86_64")
elseif (CMAKE_ASM_NASM_OBJECT_FORMAT MATCHES "elfx32")
  set(extra_flags "-m elf32_x86_64")
elseif (CMAKE_ASM_NASM_OBJECT_FORMAT MATCHES "elf(32)?")
  set(extra_flags "-m elf_i386")
else ()
  set(extra_flags "")
endif ()

set(CMAKE_ASM_NASM_CREATE_SHARED_LIBRARY
    "lld -shared ${extra_flags} <CMAKE_ASM_NASM_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>"
)
set(CMAKE_ASM_NASM_CREATE_STATIC_LIBRARY
    "lld -shared ${extra_flags} <CMAKE_ASM_NASM_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>"
)

