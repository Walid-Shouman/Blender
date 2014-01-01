/*
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * Contributor(s):
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/bmesh/tools/bmesh_data_transfer.h
 *  \ingroup bmesh
 */

#ifndef __BMESH_DATA_TRANSFER_H__
#define __BMESH_DATA_TRANSFER_H__

typedef enum ST_ShapekeyGroupMode {
	ST_REPLACE_ACTIVE_SHAPEKEY_GROUP = 1,
	ST_REPLACE_ENOUGH_SHAPEKEY_GROUPS = 2,
	ST_REPLACE_ALL_SHAPEKEY_GROUPS = 3,
	ST_APPEND_SHAPEKEY_GROUPS = 4
} ST_ShapekeyGroupMode;

typedef enum ReplaceGroupMode {
	REPLACE_ACTIVE_GROUP = 1,
	REPLACE_ENOUGH_GROUPS = 2,
	REPLACE_ALL_GROUPS = 3,
	APPEND_GROUPS = 4
} ReplaceGroupMode;

struct ReplaceLayerInfo {
	int src_lay_start;
	int src_lay_end;
	int dst_lay_start;
	int dst_lay_end;
} ReplaceLayerInfo;

typedef enum TransferMode {
	TRANSFER_BY_INDEX = 1,
	TRANSFER_BY_TOPOLOGY = 2,
	TRANSFER_BY_INTERPOLATION = 3,
} TransferMode;

bool BM_mesh_data_copy(BMesh *bm_src, BMesh* bm_dst, const struct ReplaceLayerInfo *replace_info, int type,
                        TransferMode transfer_mode, bool relative_to_target, float tmp_mat[4][4], bool use_tolerance,
                        float tolerance);

#endif /* __BMESH_DATA_TRANSFER_H__ */

