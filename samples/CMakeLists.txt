SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(sequences_example sequences_example.c)
add_dependencies(sequences_example bioc)
target_link_libraries(sequences_example bioc)

add_executable(parsing_fasta parsing_fasta.c)
add_dependencies(parsing_fasta bioc)
target_link_libraries(parsing_fasta bioc)

add_executable(bench bench.c)
add_dependencies(bench bioc)
target_link_libraries(bench bioc)
