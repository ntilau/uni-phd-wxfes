/**
 * This file has no copyright assigned and is placed in the Public Domain.
 * This file is part of the mingw-w64 runtime package.
 * No warranty is given; refer to the file DISCLAIMER.PD within this package.
 */
#ifndef _INC_NETLISTMGR
#define _INC_NETLISTMGR

#if (_WIN32_WINNT >= 0x0600)

#ifdef __cplusplus
extern "C" {
#endif

typedef enum NLM_CONNECTION_PROPERTY_CHANGE {
  NLM_CONNECTION_PROPERTY_CHANGE_AUTHENTICATION   = 0x01
} NLM_CONNECTION_PROPERTY_CHANGE;

typedef enum NLM_CONNECTIVITY {
  NLM_CONNECTIVITY_DISCONNECTED        = 0x0000,
  NLM_CONNECTIVITY_IPV4_NOTRAFFIC      = 0x0001,
  NLM_CONNECTIVITY_IPV6_NOTRAFFIC      = 0x0002,
  NLM_CONNECTIVITY_IPV4_SUBNET         = 0x0010,
  NLM_CONNECTIVITY_IPV4_LOCALNETWORK   = 0x0020,
  NLM_CONNECTIVITY_IPV4_INTERNET       = 0x0040,
  NLM_CONNECTIVITY_IPV6_SUBNET         = 0x0100,
  NLM_CONNECTIVITY_IPV6_LOCALNETWORK   = 0x0200,
  NLM_CONNECTIVITY_IPV6_INTERNET       = 0x0400
} NLM_CONNECTIVITY;

typedef enum NLM_DOMAIN_TYPE {
  NLM_DOMAIN_TYPE_NON_DOMAIN_NETWORK     = 0x0,
  NLM_DOMAIN_TYPE_DOMAIN_NETWORK         = 0x01,
  NLM_DOMAIN_TYPE_DOMAIN_AUTHENTICATED   = 0x02
} NLM_DOMAIN_TYPE;

typedef enum NLM_ENUM_NETWORK {
  NLM_ENUM_NETWORK_CONNECTED      = 0x01,
  NLM_ENUM_NETWORK_DISCONNECTED   = 0x02,
  NLM_ENUM_NETWORK_ALL            = 0x03
} NLM_ENUM_NETWORK;

typedef enum NLM_NETWORK_CATEGORY {
  NLM_NETWORK_CATEGORY_PUBLIC                 = 0x00,
  NLM_NETWORK_CATEGORY_PRIVATE                = 0x01,
  NLM_NETWORK_CATEGORY_DOMAIN_AUTHENTICATED   = 0x02
} NLM_NETWORK_CATEGORY;

typedef enum _NLM_NETWORK_CLASS {
  NLM_NETWORK_IDENTIFYING    = 0x01,
  NLM_NETWORK_IDENTIFIED     = 0x02,
  NLM_NETWORK_UNIDENTIFIED   = 0x03
} NLM_NETWORK_CLASS;

typedef enum NLM_NETWORK_PROPERTY_CHANGE {
  NLM_NETWORK_PROPERTY_CHANGE_CONNECTION       = 0x01,
  NLM_NETWORK_PROPERTY_CHANGE_DESCRIPTION      = 0x02,
  NLM_NETWORK_PROPERTY_CHANGE_NAME             = 0x04,
  NLM_NETWORK_PROPERTY_CHANGE_CATEGORY_VALUE   = 0x10
} NLM_NETWORK_PROPERTY_CHANGE;

#ifdef __cplusplus
}
#endif

#endif /*(_WIN32_WINNT >= 0x0600)*/

#endif /*_INC_NETLISTMGR*/
