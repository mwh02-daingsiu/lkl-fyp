# SPDX-License-Identifier: GPL-2.0
#
# Makefile for the linux ext2-filesystem routines.
#

obj-$(CONFIG_EXT2FYP_FS) += ext2fyp.o

ext2fyp-y := balloc.o dir.o file.o ialloc.o inode.o \
	  ioctl.o namei.o super.o symlink.o bmpt.o \
	  dup_buf.o page-io.o

ext2fyp-$(CONFIG_EXT2FYP_FS_XATTR)	 += xattr.o xattr_user.o xattr_trusted.o
ext2fyp-$(CONFIG_EXT2FYP_FS_POSIX_ACL) += acl.o
ext2fyp-$(CONFIG_EXT2FYP_FS_SECURITY)	 += xattr_security.o
