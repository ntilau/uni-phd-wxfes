/*! \file netcdf_meta.h
 *
 * Meta information for libnetcdf which can be used by other packages which
 * depend on libnetcdf.
 *
 * This file is automatically generated by the build system
 * at configure time, and contains information related to
 * how libnetcdf was built.  It will not be required to
 * include this file unless you want to probe the capabilities
 * of libnetcdf. This should ideally only happen when configuring
 * a project which depends on libnetcdf.  At configure time,
 * the dependent project can set its own macros which can be used
 * in conditionals.
 *
 * Copyright 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002,
 * 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014,
 * 2015, 2016, 2017, 2018
 * University Corporation for Atmospheric Research/Unidata.
 *
 * See \ref copyright file for more info.
 */

#ifndef NETCDF_META_H
#define NETCDF_META_H

#define NC_VERSION_MAJOR 4 /*!< netcdf-c major version. */
#define NC_VERSION_MINOR 6 /*!< netcdf-c minor version. */
#define NC_VERSION_PATCH 1 /*!< netcdf-c patch version. */
#define NC_VERSION_NOTE  "" /*!< netcdf-c note. May be blank. */

/*! netcdf-c version string.
 *
 * The format of the netcdf-c version string is as follows:
 * \code
 * NC_VERSION_MAJOR.NC_VERSION_MINOR.NC_VERSION_PATCH.NC_VERSION_NOTE
 * \endcode
 * Note that the trailing NC_VERSION_NOTE may be empty. It serves for
 * identifiers such as '-rc1', etc.
 */
#define NC_VERSION       "4.6.1"

#define NC_HAS_NC2       0 /*!< API version 2 support. */
#define NC_HAS_NC4       0 /*!< API version 4 support. */
#define NC_HAS_HDF4      0 /*!< HDF4 support. */
#define NC_HAS_HDF5      0 /*!< HDF5 support. */
#define NC_HAS_SZIP      0 /*!< szip support (HDF5 only) */
#define NC_HAS_DAP2      0 /*!< DAP2 support. */
#define NC_HAS_DAP4      0 /*!< DAP4 support. */
#define NC_HAS_BYTERANGE 
#define NC_HAS_DISKLESS  1 /*!< diskless support. */
#define NC_HAS_MMAP      0 /*!< mmap support. */
#define NC_HAS_JNA       1 /*!< jna support. */
#define NC_HAS_PNETCDF   0 /*!< PnetCDF support. */
#define NC_HAS_PARALLEL4 0 /*!< parallel IO support via HDF5 */
#define NC_HAS_PARALLEL  0 /*!< parallel IO support via HDF5 and/or PnetCDF. */

#define NC_HAS_CDF5      0  /*!< CDF5 support. */
#define NC_HAS_ERANGE_FILL  /*!< ERANGE_FILL Support */
#define NC_RELAX_COORD_BOUND 0 /*!< RELAX_COORD_BOUND */
#define NC_DISPATCH_VERSION  /*!< Dispatch table version */

#endif
