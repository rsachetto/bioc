cmake_minimum_required(VERSION 3.21)
project(bioc C)

set(CMAKE_C_STANDARD 99)

add_library(bioc src/bioc.h src/sequence.c src/sequence.h 
            src/sequence_record.c src/sequence_record.h
            src/file_utils/file_utils.c src/file_utils/file_utils.h 
            src/3d_party/hash/ht.c src/3d_party/hash/ht.h 
            src/parsing_utils.c src/parsing_utils.h)

add_subdirectory(samples)