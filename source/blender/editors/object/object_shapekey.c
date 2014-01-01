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
 * The Original Code is Copyright (C) 2001-2002 by NaN Holding BV.
 * All rights reserved.
 *
 * Contributor(s): Blender Foundation, shapekey support
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file blender/editors/object/object_shapekey.c
 *  \ingroup edobj
 */


#include <math.h>
#include <string.h>

#ifndef WIN32
#include <unistd.h>
#else
#include <io.h>
#endif   

#include "MEM_guardedalloc.h"

#include "BLI_blenlib.h"
#include "BLI_math.h"
#include "BLI_utildefines.h"

#include "DNA_curve_types.h"
#include "DNA_key_types.h"
#include "DNA_lattice_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_scene_types.h"
#include "DNA_object_types.h"

#include "BKE_context.h"
#include "BKE_depsgraph.h"
#include "BKE_key.h"
#include "BKE_library.h"
#include "BKE_main.h"
#include "BKE_object.h"
#include "BKE_lattice.h"
#include "BKE_curve.h"

#include "BLI_sys_types.h" // for intptr_t support

#include "ED_object.h"
#include "ED_mesh.h"

#include "RNA_access.h"
#include "RNA_define.h"

#include "WM_api.h"
#include "WM_types.h"

#include "object_intern.h"

#include "bmesh.h"				//todo: double check here, do we use other bmesh functios than the bmesh tools or not
#include "bmesh_tools.h"
#include "BKE_report.h"

/*********************** add shape key ***********************/

static void ED_object_shape_key_add(bContext *C, Scene *scene, Object *ob, int from_mix)
{
	KeyBlock *kb;
	if ((kb = BKE_object_insert_shape_key(scene, ob, NULL, from_mix))) {
		Key *key = BKE_key_from_object(ob);
		/* for absolute shape keys, new keys may not be added last */
		ob->shapenr = BLI_findindex(&key->block, kb) + 1;

		WM_event_add_notifier(C, NC_OBJECT | ND_DRAW, ob);
	}
}

/*********************** remove shape key ***********************/

static bool ED_object_shape_key_remove_all(Main *bmain, Object *ob)
{
	Key *key;

	key = BKE_key_from_object(ob);
	if (key == NULL)
		return false;

	switch (GS(key->from->name)) {
		case ID_ME: ((Mesh *)key->from)->key    = NULL; break;
		case ID_CU: ((Curve *)key->from)->key   = NULL; break;
		case ID_LT: ((Lattice *)key->from)->key = NULL; break;
	}

	BKE_libblock_free_us(&(bmain->key), key);

	return true;
}

static bool ED_object_shape_key_remove(Main *bmain, Object *ob)
{
	KeyBlock *kb, *rkb;
	Key *key;

	key = BKE_key_from_object(ob);
	if (key == NULL)
		return false;

	kb = BLI_findlink(&key->block, ob->shapenr - 1);

	if (kb) {
		for (rkb = key->block.first; rkb; rkb = rkb->next)
			if (rkb->relative == ob->shapenr - 1)
				rkb->relative = 0;

		BLI_remlink(&key->block, kb);
		key->totkey--;
		if (key->refkey == kb) {
			key->refkey = key->block.first;

			if (key->refkey) {
				/* apply new basis key on original data */
				switch (ob->type) {
					case OB_MESH:
						BKE_key_convert_to_mesh(key->refkey, ob->data);
						break;
					case OB_CURVE:
					case OB_SURF:
						BKE_key_convert_to_curve(key->refkey, ob->data, BKE_curve_nurbs_get(ob->data));
						break;
					case OB_LATTICE:
						BKE_key_convert_to_lattice(key->refkey, ob->data);
						break;
				}
			}
		}
			
		if (kb->data) MEM_freeN(kb->data);
		MEM_freeN(kb);

		if (ob->shapenr > 1) {
			ob->shapenr--;
		}
	}
	
	if (key->totkey == 0) {
		switch (GS(key->from->name)) {
			case ID_ME: ((Mesh *)key->from)->key    = NULL; break;
			case ID_CU: ((Curve *)key->from)->key   = NULL; break;
			case ID_LT: ((Lattice *)key->from)->key = NULL; break;
		}

		BKE_libblock_free_us(&(bmain->key), key);
	}

	return true;
}

static bool object_shape_key_mirror(bContext *C, Object *ob,
                                    int *r_totmirr, int *r_totfail, bool use_topology)
{
	KeyBlock *kb;
	Key *key;
	int totmirr = 0, totfail = 0;

	*r_totmirr = *r_totfail = 0;

	key = BKE_key_from_object(ob);
	if (key == NULL)
		return 0;
	
	kb = BLI_findlink(&key->block, ob->shapenr - 1);

	if (kb) {
		char *tag_elem = MEM_callocN(sizeof(char) * kb->totelem, "shape_key_mirror");


		if (ob->type == OB_MESH) {
			Mesh *me = ob->data;
			MVert *mv;
			int i1, i2;
			float *fp1, *fp2;
			float tvec[3];

			mesh_octree_table(ob, NULL, NULL, 's');

			for (i1 = 0, mv = me->mvert; i1 < me->totvert; i1++, mv++) {
				i2 = mesh_get_x_mirror_vert(ob, i1, use_topology);
				if (i2 == i1) {
					fp1 = ((float *)kb->data) + i1 * 3;
					fp1[0] = -fp1[0];
					tag_elem[i1] = 1;
					totmirr++;
				}
				else if (i2 != -1) {
					if (tag_elem[i1] == 0 && tag_elem[i2] == 0) {
						fp1 = ((float *)kb->data) + i1 * 3;
						fp2 = ((float *)kb->data) + i2 * 3;

						copy_v3_v3(tvec,    fp1);
						copy_v3_v3(fp1, fp2);
						copy_v3_v3(fp2, tvec);

						/* flip x axis */
						fp1[0] = -fp1[0];
						fp2[0] = -fp2[0];
						totmirr++;
					}
					tag_elem[i1] = tag_elem[i2] = 1;
				}
				else {
					totfail++;
				}
			}

			mesh_octree_table(ob, NULL, NULL, 'e');
		}
		else if (ob->type == OB_LATTICE) {
			Lattice *lt = ob->data;
			int i1, i2;
			float *fp1, *fp2;
			int u, v, w;
			/* half but found up odd value */
			const int pntsu_half = (lt->pntsu / 2) + (lt->pntsu % 2);

			/* currently editmode isn't supported by mesh so
			 * ignore here for now too */

			/* if (lt->editlatt) lt = lt->editlatt->latt; */

			for (w = 0; w < lt->pntsw; w++) {
				for (v = 0; v < lt->pntsv; v++) {
					for (u = 0; u < pntsu_half; u++) {
						int u_inv = (lt->pntsu - 1) - u;
						float tvec[3];
						if (u == u_inv) {
							i1 = BKE_lattice_index_from_uvw(lt, u, v, w);
							fp1 = ((float *)kb->data) + i1 * 3;
							fp1[0] = -fp1[0];
							totmirr++;
						}
						else {
							i1 = BKE_lattice_index_from_uvw(lt, u, v, w);
							i2 = BKE_lattice_index_from_uvw(lt, u_inv, v, w);

							fp1 = ((float *)kb->data) + i1 * 3;
							fp2 = ((float *)kb->data) + i2 * 3;

							copy_v3_v3(tvec, fp1);
							copy_v3_v3(fp1, fp2);
							copy_v3_v3(fp2, tvec);
							fp1[0] = -fp1[0];
							fp2[0] = -fp2[0];
							totmirr++;
						}
					}
				}
			}
		}

		MEM_freeN(tag_elem);
	}
	
	*r_totmirr = totmirr;
	*r_totfail = totfail;

	DAG_id_tag_update(&ob->id, OB_RECALC_DATA);
	WM_event_add_notifier(C, NC_OBJECT | ND_DRAW, ob);

	return 1;
}

/*********************** transfer shape key ***********************/

static EnumPropertyItem replace_mode_item[] = {
    {REPLACE_ACTIVE_GROUP,
	 "REPLACE_ACTIVE_GROUP", 0, "Active", "Overwrite active group only"},
    {REPLACE_ENOUGH_GROUPS,
	 "REPLACE_ENOUGH_GROUPS", 0, "Enough", "Overwrite source groups only as needed"},
    {REPLACE_ALL_GROUPS,
	 "REPLACE_ALL_GROUPS", 0, "All", "Overwrite all groups"},
    {APPEND_GROUPS,
	 "APPEND_GROUPS", 0, "Append", "Add groups without overwriting"},
	{0, NULL, 0, NULL, NULL}
};

typedef enum ST_FromToActive {
	ST_FROM_ACTIVE = 1,
	ST_TO_ACTIVE = 2
} ST_FromToActive;

static EnumPropertyItem ST_from_to_active[] = {
    {ST_FROM_ACTIVE,
	 "ST_FROM_ACTIVE", 0, "From active", "Transfer to different objects"},
    {ST_TO_ACTIVE,
	 "ST_TO_ACTIVE", 0, "To active", "Better to faster tweek the output"},
    {0, NULL, 0, NULL, NULL}
};

static EnumPropertyItem transfer_mode_item[] = {
    {TRANSFER_BY_INDEX,
	 "TRANSFER_BY_INDEX", 0, "By index", "copy between identical indices meshes"},
    {TRANSFER_BY_TOPOLOGY,
	 "TRANSFER_BY_TOPOLOGY", 0, "By topology", "use if the same topology with different indices"},
    {TRANSFER_BY_INTERPOLATION,
	 "TRANSFER_BY_INTERPOLATION", 0, "By interpolation", "interpolate for different topologies"},
	{0, NULL, 0, NULL, NULL}
};

static bool ED_object_shape_key_transfer(Object *ob_dst, Object *ob_src, bContext *C, Scene *scene, wmOperator *op)
{
	Mesh *me_dst, *me_src;
	BMesh *bm_dst, *bm_src;

	int *act_shapekey_lay = MEM_mallocN(sizeof(*act_shapekey_lay) * 2, "act_shapekey_lay object_shapekey.c");

	bool relative_to_target = RNA_boolean_get(op->ptr, "rel_to_target");
	ReplaceGroupMode replace_mode = RNA_enum_get(op->ptr, "replace_mode");
	bool use_tolerance = RNA_boolean_get(op->ptr, "use_tol");
	float tolerance2 = RNA_float_get(op->ptr, "tol");
	PropertyRNA *tolerance_prop = RNA_struct_find_property(op->ptr, "tol");
	TransferMode transfer_mode = RNA_enum_get(op->ptr, "transfer_mode");

	struct ReplaceLayerInfo replace_info;

	int i = 0;

	Main *bmain = CTX_data_main(C);

	KeyBlock *kb_src, *kb_dst;

	//----- definitions for the raycasting
	float tmp_mat[4][4];
	//===== end of raycasting definitions

	int active_dst, active_src;
	int num_src_lay, num_dst_lay;

	//----- raycasting assignments
	//Prepare transformation matrix.
	invert_m4_m4(ob_src->imat, ob_src->obmat);
	mul_m4_m4m4(tmp_mat, ob_src->imat, ob_dst->obmat);
	//===== end of RC assignments

	if(use_tolerance == false)
		RNA_def_property_flag(tolerance_prop, PROP_HIDDEN);
	else
		RNA_def_property_clear_flag(tolerance_prop, PROP_HIDDEN);

	me_dst = ob_dst->data;
	me_src = ob_src->data;

	//allocate space
	bm_src = BM_mesh_create(&bm_mesh_allocsize_default);
	bm_dst = BM_mesh_create(&bm_mesh_allocsize_default);

	if (me_src->key != NULL) {
		if (me_src->key->totkey < 2) {
			//the source should have at least a basis and one more layer
			BKE_report(op->reports, RPT_ERROR,
			           "Transfer failed (source mesh should have -at least- the basis and another layer)");
			return false;
		}
	}

	else {
		BKE_report(op->reports, RPT_ERROR,
		           "Transfer failed no shapekeys were found (source mesh should have -at least- the basis and another layer)");
		return false;
	}

	num_src_lay = me_src->key->totkey;
	num_dst_lay = (me_dst->key == NULL) ? 0: me_dst->key->totkey;

	//we'll also tell the copy function to start copying from the first shapekey after basis
	if (replace_mode == REPLACE_ENOUGH_GROUPS) {
		//add the dst basis if not found
		if (me_dst->key == NULL) {
			ED_object_shape_key_add(C, scene, ob_dst, false);
		}

		kb_src = me_src->key->block.first;
		kb_dst = me_dst->key->block.first;

		//add layers as needed
		while (me_dst->key->totkey < me_src->key->totkey) {
			ED_object_shape_key_add(C, scene, ob_dst, false);
		}

		//copy the names
		for (i = 1; i < me_src->key->totkey; ++i) {
			kb_src = kb_src->next;
			kb_dst = kb_dst->next;

			BLI_strncpy(kb_dst->name, kb_src->name, sizeof(kb_src->name));
		}

		replace_info.src_lay_start = 1;
		replace_info.src_lay_end = num_src_lay - 1;
		replace_info.dst_lay_start = 1;
		replace_info.dst_lay_end = num_src_lay - 1;
	}

	//we'll tell the copy function to start copying from # of source keys from the end of the dst keys
	else if (replace_mode == APPEND_GROUPS)
	{
		//add the dst basis if not found
		if (me_dst->key == NULL) {
			ED_object_shape_key_add(C, scene, ob_dst, false);
			num_dst_lay++;	//we just can't let the basis be copied into
		}

		kb_src = me_src->key->block.first;
		kb_dst = me_dst->key->block.last;

		//skip the src basis while appending
		for (i = 1; i < me_src->key->totkey; ++i) {
			kb_src = kb_src->next;

			ED_object_shape_key_add(C, scene, ob_dst, false);
			kb_dst = kb_dst->next;

			//rename each appended layer
			BLI_strncpy(kb_dst->name, kb_src->name, sizeof(kb_src->name));
		}

		replace_info.src_lay_start = 1;
		replace_info.src_lay_end = num_src_lay - 1;
		replace_info.dst_lay_start = num_dst_lay;
		replace_info.dst_lay_end = num_dst_lay + num_src_lay - 2;
	}
	//same message to be sent to the copy function as ST_REPLACE_ENOUGH_VERTEX_GROUPS
	else if (replace_mode == REPLACE_ALL_GROUPS)
	{
		//add the dst basis if not found

		if (me_dst->key == NULL) {
			ED_object_shape_key_add(C, scene, ob_dst, false);
		}

		kb_src = me_src->key->block.first;
		kb_dst = me_dst->key->block.first;

		//start from 1 to skip deleting the basis
		while (me_dst->key->totkey > 1) {
			ED_object_shape_key_remove(bmain , ob_dst);
		}

		//add layers as needed
		while (me_dst->key->totkey < me_src->key->totkey) {
			ED_object_shape_key_add(C, scene, ob_dst, false);
		}

		//copy the names
		for (i = 1; i < me_src->key->totkey; ++i) {
			kb_src = kb_src->next;
			kb_dst = kb_dst->next;

			BLI_strncpy(kb_dst->name, kb_src->name, sizeof(kb_src->name));
		}

		replace_info.src_lay_start = 1;
		replace_info.src_lay_end = num_src_lay - 1;
		replace_info.dst_lay_start = 1;
		replace_info.dst_lay_end = num_src_lay - 1;
	}

	else if (replace_mode == REPLACE_ACTIVE_GROUP) {

		active_src = ob_src->shapenr;
		active_dst = ob_dst->shapenr;

		if (active_src > 1) {

			//find the source
			kb_src = BLI_findlink(&me_src->key->block, active_src - 1);

			if (active_dst == 0) {	//empty destination
				ED_object_shape_key_add(C, scene, ob_dst, false);
				active_dst++;
			}

			if (active_dst == 1) { //destination's basis is selected
				ED_object_shape_key_add(C, scene, ob_dst, false);
				kb_dst = ((KeyBlock *) me_dst->key->block.first)->next;

				BLI_strncpy(kb_dst->name, kb_src->name, sizeof(kb_src->name));
				active_dst++;
			}

			else {		//another layer is selected
				kb_dst = BLI_findlink(&me_dst->key->block, active_dst - 1);

				BLI_strncpy(kb_dst->name, kb_src->name, sizeof(kb_src->name));
			}

			replace_info.src_lay_start = active_src - 1;
			replace_info.src_lay_end = replace_info.src_lay_start;
			replace_info.dst_lay_start = active_dst - 1;	//fixing the indices
			replace_info.dst_lay_end = replace_info.dst_lay_start;
		}

		else {
			BKE_report(op->reports, RPT_ERROR,
			           "Transfer failed (The active shapekey group isn't a valid one ensure it's not the basis)");
			return false;
		}
	}

	BM_mesh_bm_from_me(bm_src, me_src, TRUE, true, 0);	//TRUE -> should transfer shapekeys too!!
	BM_mesh_bm_from_me(bm_dst, me_dst, TRUE, true, 0);

	if (!BM_mesh_data_copy(bm_src, bm_dst, &replace_info, CD_SHAPEKEY, transfer_mode, relative_to_target, tmp_mat,
	                       use_tolerance, tolerance2)) {
	return false;
	}

	//transfer the BMesh back to Mesh
	BM_mesh_bm_to_me(bm_src, me_src, FALSE);
	BM_mesh_bm_to_me(bm_dst, me_dst, TRUE);

	//free the BMesh
	BM_mesh_free(bm_src);
	BM_mesh_free(bm_dst);
	return true;
}

/********************** shape key operators *********************/

static int shape_key_mode_poll(bContext *C)
{
	Object *ob = ED_object_context(C);
	ID *data = (ob) ? ob->data : NULL;
	return (ob && !ob->id.lib && data && !data->lib && ob->mode != OB_MODE_EDIT);
}

static int shape_key_mode_exists_poll(bContext *C)
{
	Object *ob = ED_object_context(C);
	ID *data = (ob) ? ob->data : NULL;

	/* same as shape_key_mode_poll */
	return (ob && !ob->id.lib && data && !data->lib && ob->mode != OB_MODE_EDIT) &&
	       /* check a keyblock exists */
	       (BKE_keyblock_from_object(ob) != NULL);
}

static int shape_key_poll(bContext *C)
{
	Object *ob = ED_object_context(C);
	ID *data = (ob) ? ob->data : NULL;
	return (ob && !ob->id.lib && data && !data->lib);
}

static int shape_key_add_exec(bContext *C, wmOperator *op)
{
	Scene *scene = CTX_data_scene(C);
	Object *ob = ED_object_context(C);
	int from_mix = RNA_boolean_get(op->ptr, "from_mix");

	ED_object_shape_key_add(C, scene, ob, from_mix);

	return OPERATOR_FINISHED;
}

void OBJECT_OT_shape_key_add(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "Add Shape Key";
	ot->idname = "OBJECT_OT_shape_key_add";
	ot->description = "Add shape key to the object";
	
	/* api callbacks */
	ot->poll = shape_key_mode_poll;
	ot->exec = shape_key_add_exec;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

	/* properties */
	RNA_def_boolean(ot->srna, "from_mix", 1, "From Mix", "Create the new shape key from the existing mix of keys");
}

static int shape_key_remove_exec(bContext *C, wmOperator *op)
{
	Main *bmain = CTX_data_main(C);
	Object *ob = ED_object_context(C);
	bool changed = false;

	if (RNA_boolean_get(op->ptr, "all")) {
		changed = ED_object_shape_key_remove_all(bmain, ob);
	}
	else {
		changed = ED_object_shape_key_remove(bmain, ob);
	}

	if (changed) {
		DAG_id_tag_update(&ob->id, OB_RECALC_DATA);
		WM_event_add_notifier(C, NC_OBJECT | ND_DRAW, ob);

		return OPERATOR_FINISHED;
	}
	else {
		return OPERATOR_CANCELLED;
	}
}

void OBJECT_OT_shape_key_remove(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "Remove Shape Key";
	ot->idname = "OBJECT_OT_shape_key_remove";
	ot->description = "Remove shape key from the object";
	
	/* api callbacks */
	ot->poll = shape_key_mode_poll;
	ot->poll = shape_key_mode_exists_poll;
	ot->exec = shape_key_remove_exec;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

	/* properties */
	RNA_def_boolean(ot->srna, "all", 0, "All", "Remove all shape keys");
}

static int shape_key_clear_exec(bContext *C, wmOperator *UNUSED(op))
{
	Object *ob = ED_object_context(C);
	Key *key = BKE_key_from_object(ob);
	KeyBlock *kb = BKE_keyblock_from_object(ob);

	if (!key || !kb)
		return OPERATOR_CANCELLED;
	
	for (kb = key->block.first; kb; kb = kb->next)
		kb->curval = 0.0f;

	DAG_id_tag_update(&ob->id, OB_RECALC_DATA);
	WM_event_add_notifier(C, NC_OBJECT | ND_DRAW, ob);
	
	return OPERATOR_FINISHED;
}

void OBJECT_OT_shape_key_clear(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "Clear Shape Keys";
	ot->description = "Clear weights for all shape keys";
	ot->idname = "OBJECT_OT_shape_key_clear";
	
	/* api callbacks */
	ot->poll = shape_key_poll;
	ot->exec = shape_key_clear_exec;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

/* starting point and step size could be optional */
static int shape_key_retime_exec(bContext *C, wmOperator *UNUSED(op))
{
	Object *ob = ED_object_context(C);
	Key *key = BKE_key_from_object(ob);
	KeyBlock *kb = BKE_keyblock_from_object(ob);
	float cfra = 0.0f;

	if (!key || !kb)
		return OPERATOR_CANCELLED;

	for (kb = key->block.first; kb; kb = kb->next)
		kb->pos = (cfra += 0.1f);

	DAG_id_tag_update(&ob->id, OB_RECALC_DATA);
	WM_event_add_notifier(C, NC_OBJECT | ND_DRAW, ob);

	return OPERATOR_FINISHED;
}

void OBJECT_OT_shape_key_retime(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "Re-Time Shape Keys";
	ot->description = "Resets the timing for absolute shape keys";
	ot->idname = "OBJECT_OT_shape_key_retime";

	/* api callbacks */
	ot->poll = shape_key_poll;
	ot->exec = shape_key_retime_exec;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;
}

static int shape_key_mirror_exec(bContext *C, wmOperator *op)
{
	Object *ob = ED_object_context(C);
	int totmirr = 0, totfail = 0;
	bool use_topology = RNA_boolean_get(op->ptr, "use_topology");

	if (!object_shape_key_mirror(C, ob, &totmirr, &totfail, use_topology))
		return OPERATOR_CANCELLED;

	ED_mesh_report_mirror(op, totmirr, totfail);

	return OPERATOR_FINISHED;
}

void OBJECT_OT_shape_key_mirror(wmOperatorType *ot)
{
	/* identifiers */
	ot->name = "Mirror Shape Key";
	ot->idname = "OBJECT_OT_shape_key_mirror";
	ot->description = "Mirror the current shape key along the local X axis";

	/* api callbacks */
	ot->poll = shape_key_mode_poll;
	ot->exec = shape_key_mirror_exec;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

	/* properties */
	RNA_def_boolean(ot->srna, "use_topology", 0, "Topology Mirror",
	                "Use topology based mirroring (for when both sides of mesh have matching, unique topology)");
}


static int shape_key_move_exec(bContext *C, wmOperator *op)
{
	Object *ob = ED_object_context(C);

	int type = RNA_enum_get(op->ptr, "type");
	Key *key = BKE_key_from_object(ob);

	if (key) {
		KeyBlock *kb, *kb_other;
		int shapenr_act = ob->shapenr - 1;
		int shapenr_swap = shapenr_act + type;
		kb = BLI_findlink(&key->block, shapenr_act);

		if ((type == -1 && kb->prev == NULL) || (type == 1 && kb->next == NULL)) {
			return OPERATOR_CANCELLED;
		}

		for (kb_other = key->block.first; kb_other; kb_other = kb_other->next) {
			if (kb_other->relative == shapenr_act) {
				kb_other->relative += type;
			}
			else if (kb_other->relative == shapenr_swap) {
				kb_other->relative -= type;
			}
		}

		if (type == -1) {
			/* move back */
			kb_other = kb->prev;
			BLI_remlink(&key->block, kb);
			BLI_insertlinkbefore(&key->block, kb_other, kb);
			ob->shapenr--;
		}
		else {
			/* move next */
			kb_other = kb->next;
			BLI_remlink(&key->block, kb);
			BLI_insertlinkafter(&key->block, kb_other, kb);
			ob->shapenr++;
		}

		SWAP(float, kb_other->pos, kb->pos); /* for absolute shape keys */

		/* First key is refkey, matches interface and BKE_key_sort */
		key->refkey = key->block.first;
	}

	DAG_id_tag_update(&ob->id, OB_RECALC_DATA);
	WM_event_add_notifier(C, NC_OBJECT | ND_DRAW, ob);

	return OPERATOR_FINISHED;
}

void OBJECT_OT_shape_key_move(wmOperatorType *ot)
{
	static EnumPropertyItem slot_move[] = {
		{-1, "UP", 0, "Up", ""},
		{1, "DOWN", 0, "Down", ""},
		{0, NULL, 0, NULL, NULL}
	};

	/* identifiers */
	ot->name = "Move Shape Key";
	ot->idname = "OBJECT_OT_shape_key_move";
	ot->description = "Move the active shape key up/down in the list";

	/* api callbacks */
	ot->poll = shape_key_mode_poll;
	ot->exec = shape_key_move_exec;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;

	RNA_def_enum(ot->srna, "type", slot_move, 0, "Type", "");
}

static int shape_key_transfer_exec(bContext *C, wmOperator *op)
{
	Scene *scene = CTX_data_scene(C);
	Object *ob_act = CTX_data_active_object(C);
	int fail = 0;

	bool transfer_first_to_act = true;

	ST_FromToActive from_active = RNA_enum_get(op->ptr, "from_to_active");

	/* Macro to loop through selected objects.*/
	CTX_DATA_BEGIN (C, Object *, ob_slc, selected_editable_objects)
	{
		//if the selected isn't the active object
		if (ob_act != ob_slc) {

			if (from_active == ST_TO_ACTIVE) {

				//if many objects were selected within this mode ... we should copy only from the first
				//notice that ob_slc priority isn't set by order of selection!
				if (transfer_first_to_act == true) {
					transfer_first_to_act = false;

					if (!ED_object_shape_key_transfer(ob_act, ob_slc, C, scene, op)) {
						fail++;
						}
				}

			}
			else {		//copy from the active to all the other selected
				if (!ED_object_shape_key_transfer(ob_slc, ob_act, C, scene, op)) {
					fail++;
				}
			}
		}
	}

	////ported from transfer weights
	/* Event notifiers for correct display of data.*/
	DAG_id_tag_update(&ob_slc->id, OB_RECALC_DATA);
	WM_event_add_notifier(C, NC_OBJECT | ND_DRAW, ob_slc);
	WM_event_add_notifier(C, NC_GEOM | ND_DATA, ob_slc->data);

	CTX_DATA_END;

	if (fail != 0) {
		return OPERATOR_CANCELLED;
	}
	else {
		return OPERATOR_FINISHED;
	}
}

void OBJECT_OT_shape_key_transfer_new(wmOperatorType *ot)
{

	/* identifiers */
	ot->name = "Transfer Shape Key (new)";
	ot->idname = "OBJECT_OT_shape_key_transfer_new";
	ot->description = "Transfer shapekey groups to the selected objects";

	/* api callbacks */
	ot->poll = shape_key_poll;			//don't know how to edit this yet!!
	ot->exec = shape_key_transfer_exec;

	/* flags */
	ot->flag = OPTYPE_REGISTER | OPTYPE_UNDO;	//not revised!!

	/* properties */
	RNA_def_boolean(ot->srna, "rel_to_target", false,
	                "Relative to target", "select this if you want the transfer to be relative to the target");
	RNA_def_enum(ot->srna, "replace_mode", replace_mode_item, 2,
	             "Replace/Append", "define which groups to move");
	RNA_def_enum(ot->srna, "from_to_active", ST_from_to_active, 2, "From/To active object",
	             "Choose the transfer direction");
	RNA_def_boolean(ot->srna, "use_tol", false, "Use Tolerance",
	                "use a tolerance less than infinity to search for the nearest source faces");
	RNA_def_float(ot->srna, "tol", 1, 0, FLT_MAX, "Tolerance",
	              "Overwrite the search area to be a value other than infinity; useful for partial transfer", 0, 1000);
	RNA_def_enum(ot->srna, "transfer_mode", transfer_mode_item, 1,
	             "index, topology or interpolate", "define which groups to move");
}

