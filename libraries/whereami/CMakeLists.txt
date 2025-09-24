##########################################
# whereami
##########################################

add_library(
  whereami
  whereami.c 
  )

target_include_directories(whereami PUBLIC ./)