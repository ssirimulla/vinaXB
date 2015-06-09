/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#include "naive_non_cache.h"
#include "curl.h"

naive_non_cache::naive_non_cache(const precalculate* p_) : p(p_) {}
fl naive_non_cache::eval(const model& m, fl v) const { // needs m.coords
	fl e = 0;
	const fl cutoff_sqr = p->cutoff_sqr();

	sz n = num_atom_types(p->atom_typing_used());

	VINA_FOR(i, m.num_movable_atoms()) {
		fl this_e = 0;
		const atom& a = m.atoms[i];
		sz t1 = a.get(p->atom_typing_used());
		if(t1 >= n) continue;
		const vec& a_coords = m.coords[i];

		VINA_FOR_IN(j, m.grid_atoms) {
			const atom& b = m.grid_atoms[j];
			sz t2 = b.get(p->atom_typing_used());
			if(t2 >= n) continue;
			vec r_ba; r_ba = a_coords - b.coords;
			fl r2 = sqr(r_ba);
            fl theta;
            if (xs_hal_any_bond_possible(a.xs, b.xs)) {     //if there is a halogen bond possible, we need the angle
                if (xs_is_halogen(a.xs)) {                  //a is the halogen, b is the acceptor
                    const std::vector<bond>& bonds = a.bonds;
                    VINA_FOR_IN(j, bonds) {
                        const bond& bs = bonds[j];
                        const atom& c = m.get_atom(bs.connected_atom_index);
                        if(c.el == EL_TYPE_C){              //the halogen needs to be connected to a carbon
                            vec C_coords = m.coords[bs.connected_atom_index.i];     //coords of carbon
                            vec v1 = a_coords - C_coords;  //vector between halogen and carbon
                            vec v2 = a_coords - b.coords;  //vector between halogen and acceptor
                            theta = angle(v1, v2);
                        }
                    }
                }
                else {                                  //b is the halogen, a is the acceptor (might never happen?)
                    const std::vector<bond>& bonds = b.bonds;
                    VINA_FOR_IN(j, bonds) {
                        const bond& bs = bonds[j];
                        const atom& c = m.get_atom(bs.connected_atom_index);
                        if(c.el == EL_TYPE_C){              //the halogen needs to be connected to a carbon
                            vec v1 = b.coords - c.coords;  //vector between halogen and carbon
                            vec v2 = b.coords - a_coords;  //vector between halogen and acceptor
                            fl theta = angle(v1, v2);
                        }
                    }
                }
            }
            else {
                theta = 180;
            }
			if(r2 < cutoff_sqr) {
				sz type_pair_index = get_type_pair_index(p->atom_typing_used(), a, b);
                this_e +=  p->eval_fast(type_pair_index, r2, theta);
			}
		}
		curl(this_e, v);
		e += this_e;
	}
	return e;
}
