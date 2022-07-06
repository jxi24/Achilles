add_library(sherpa SHARED
    SherpaMEs.cc
    Channels3.cc
    Channels1.cc
)
target_include_directories(sherpa PUBLIC ${SHERPA_INCLUDE_DIRS})
target_link_libraries(sherpa PRIVATE project_options
                             PUBLIC ${SHERPA_LIBRARIES} utilities -L/Users/sherry/Desktop/Fermilab/SHERPA-MC-2.3.0/lib/SHERPA-MC/ -lModelMain -lToolsOrg -lPhasicChannels -lToolsPhys -lPhasicMain -lComix -lPhasicProcess -lPhasicSelectors spdlog::spdlog yaml::cpp)

if(ENABLE_GZIP)
    target_link_libraries(sherpa PUBLIC ZLIB::ZLIB)
endif()
