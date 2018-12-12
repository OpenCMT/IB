# Target for RPMs creation

set(PKGREL 1)
set(ITEMS_FOR_SOURCE AUTHORS
                     CMake
                     CMakeConfig.in.h
                     CMakeLists.txt
                     IB.config
                     IBConfig.cmake.in
                     IBConfigVersion.cmake.in
                     Jenkinsfile
                     LICENSE
                     recipes
                     src
                     utils)

add_custom_target(rpm 
                  COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/rpmbuild/BUILD
                                   ${CMAKE_CURRENT_BINARY_DIR}/rpmbuild/RPMS
                                   ${CMAKE_CURRENT_BINARY_DIR}/rpmbuild/SOURCES
                                   ${CMAKE_CURRENT_BINARY_DIR}/rpmbuild/SPECS
                                   ${CMAKE_CURRENT_BINARY_DIR}/rpmbuild/SRPMS
                  COMMAND tar -zcf ${CMAKE_CURRENT_BINARY_DIR}/rpmbuild/SOURCES/cmt-ib-${PACKAGE_VERSION}.tar.gz ${ITEMS_FOR_SOURCE}
                  COMMAND sed -e 's|@PKGVERSION@|${PACKAGE_VERSION}|g'
                              -e 's|@PKGRELEASE@|${PKGREL}|g' 
                              CMake/cmt-ib.spec.in > ${CMAKE_CURRENT_BINARY_DIR}/rpmbuild/SPECS/cmt-ib.spec
                  COMMAND QA_SKIP_BUILD_ROOT=1 rpmbuild -ba --define '_topdir ${CMAKE_CURRENT_BINARY_DIR}/rpmbuild'
                          ${CMAKE_CURRENT_BINARY_DIR}/rpmbuild/SPECS/cmt-ib.spec
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

