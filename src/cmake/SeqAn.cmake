if (SEQAN_ROOT)
	if (NOT EXISTS ${SEQAN_ROOT}/share/cmake/Modules/FindSeqAn.cmake)
		message ( FATAL_ERROR "SEQAN_ROOT was specified but '${SEQAN_ROOT}/share/cmake/Modules/FindSeqAn.cmake' does not exist." )
	endif()
else()
	set ( SEQAN_URL "https://github.com/skrakau/seqan/releases/download/seqan-v2.2.0_mod/seqan-library-2.2.0_mod.zip")
	set ( SEQAN_MD5 "7758a6b8f3000e07d580cfdd4cc57f2f")
	set ( SEQAN_ZIP_OUT ${CMAKE_CURRENT_BINARY_DIR}/seqan-library-2.2.0_mod.zip )
	set ( SEQAN_ROOT ${CMAKE_CURRENT_BINARY_DIR}/seqan-library-2.2.0_mod )

	if (NOT EXISTS ${SEQAN_ROOT}/share/cmake/Modules/FindSeqAn.cmake )
		# Download zip file
		message ("Downloading ${SEQAN_URL}")
		file (DOWNLOAD ${SEQAN_URL} ${SEQAN_ZIP_OUT}
			EXPECTED_MD5 ${SEQAN_MD5}
			SHOW_PROGRESS STATUS status)
		list ( GET status 0 ret )
		list ( GET status 0 str)
		if ( NOT ret EQUAL 0)
			message (FATAL_ERROR "Download failed")
		endif()
		# Unpack zip file
		message ("Unpacking ${SEQAN_ZIP_OUT}")
        execute_process(COMMAND cmake -E tar zxf ${SEQAN_ZIP_OUT} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
		# Remove zip file
		if ( EXISTS ${SEQAN_ZIP_OUT} )
			file (REMOVE ${SEQAN_ZIP_OUT})
		endif()
	endif()
	# Check if FindSeqAn.cmake can be found wher it should
	if (NOT EXISTS ${SEQAN_ROOT}/share/cmake/Modules/FindSeqAn.cmake )
		message (FATAL_ERROR "Failed to download and unpack '${SEQAN_URL}'")
	endif()
endif ()

LIST ( APPEND CMAKE_MODULE_PATH ${SEQAN_ROOT}/share/cmake/Modules/ )
SET( SEQAN_INCLUDE_PATH ${SEQAN_ROOT}/include/ )

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB REQUIRED)

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)
