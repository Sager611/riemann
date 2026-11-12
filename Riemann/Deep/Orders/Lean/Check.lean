import Riemann.Deep.Orders.Lean.Main

/-! Axiom check: every main theorem uses only propext, Classical.choice, Quot.sound. -/

#print axioms Orders.check_iff
#print axioms Orders.perfect_extend
#print axioms Orders.bulk_steps
#print axioms Orders.bulk_vals
#print axioms Orders.patch_sound
#print axioms Orders.works_of
#print axioms Orders.reach_sound
#print axioms Orders.noPerfect_sound
#print axioms Orders.bound_of_reach
#print axioms Orders.corollary
#print axioms Orders.D_pxmd.works
#print axioms Orders.D_pxdm.works
#print axioms Orders.D_pdxm.works
#print axioms Orders.D_pdmx.works
#print axioms Orders.D_xpmd.works
#print axioms Orders.D_xpdm.works
#print axioms Orders.D_xdpm.works
#print axioms Orders.D_xdmp.works
#print axioms Orders.D_dxpm.works
#print axioms Orders.D_dpxm.works
#print axioms Orders.D_dpmx.works
#print axioms Orders.fail_pmxd
#print axioms Orders.fail_pmdx
#print axioms Orders.fail_xmpd
#print axioms Orders.fail_xmdp
#print axioms Orders.fail_mpxd
#print axioms Orders.fail_mpdx
#print axioms Orders.fail_mxpd
#print axioms Orders.fail_mxdp
#print axioms Orders.fail_mdpx
#print axioms Orders.fail_mdxp
#print axioms Orders.fail_dxmp
#print axioms Orders.fail_dmpx
#print axioms Orders.fail_dmxp
#print axioms Orders.small_pxmd
#print axioms Orders.small_pxdm
#print axioms Orders.small_pmxd
#print axioms Orders.small_pmdx
#print axioms Orders.small_pdxm
#print axioms Orders.small_pdmx
#print axioms Orders.small_xpmd
#print axioms Orders.small_xpdm
#print axioms Orders.small_xmpd
#print axioms Orders.small_xmdp
#print axioms Orders.small_xdpm
#print axioms Orders.small_xdmp
#print axioms Orders.small_mpxd
#print axioms Orders.small_mpdx
#print axioms Orders.small_mxpd
#print axioms Orders.small_mxdp
#print axioms Orders.small_mdpx
#print axioms Orders.small_mdxp
#print axioms Orders.small_dpxm
#print axioms Orders.small_dpmx
#print axioms Orders.small_dxpm
#print axioms Orders.small_dxmp
#print axioms Orders.small_dmpx
#print axioms Orders.small_dmxp
#print axioms Orders.thm1_pxmd
#print axioms Orders.thm1_pxdm
#print axioms Orders.thm1_pmxd
#print axioms Orders.thm1_pmdx
#print axioms Orders.thm1_pdxm
#print axioms Orders.thm1_pdmx
#print axioms Orders.thm1_xpmd
#print axioms Orders.thm1_xpdm
#print axioms Orders.thm1_xmpd
#print axioms Orders.thm1_xmdp
#print axioms Orders.thm1_xdpm
#print axioms Orders.thm1_xdmp
#print axioms Orders.thm1_mpxd
#print axioms Orders.thm1_mpdx
#print axioms Orders.thm1_mxpd
#print axioms Orders.thm1_mxdp
#print axioms Orders.thm1_mdpx
#print axioms Orders.thm1_mdxp
#print axioms Orders.thm1_dpxm
#print axioms Orders.thm1_dpmx
#print axioms Orders.thm1_dxpm
#print axioms Orders.thm1_dxmp
#print axioms Orders.thm1_dmpx
#print axioms Orders.thm1_dmxp

/-! Sanity checks of the checker on the paper's Table (v2 l. 110–119, perfect for N = 9). -/
example : Orders.check ⟨.add, .mul, .sub, .div⟩ 9
    [⟨2, 2, 4⟩, ⟨2, 4, 8⟩, ⟨8, 2, 6⟩, ⟨6, 2, 3⟩, ⟨2, 3, 5⟩, ⟨3, 3, 9⟩, ⟨9, 2, 7⟩, ⟨2, 2, 1⟩] = true := by
  decide
/-- operand order of − matters: 2 − 8 is not 6 -/
example : Orders.check ⟨.add, .mul, .sub, .div⟩ 9
    [⟨2, 2, 4⟩, ⟨2, 4, 8⟩, ⟨2, 8, 6⟩, ⟨6, 2, 3⟩, ⟨2, 3, 5⟩, ⟨3, 3, 9⟩, ⟨9, 2, 7⟩, ⟨2, 2, 1⟩] = false := by
  decide
/-- the schedule matters: the same moves under × + − ÷ fail at move 1 -/
example : Orders.check ⟨.mul, .add, .sub, .div⟩ 9
    [⟨2, 2, 4⟩, ⟨2, 4, 8⟩, ⟨8, 2, 6⟩, ⟨6, 2, 3⟩, ⟨2, 3, 5⟩, ⟨3, 3, 9⟩, ⟨9, 2, 7⟩, ⟨2, 2, 1⟩] = false := by
  decide
/-- one move short: 7 moves cannot be perfect for N = 9 -/
example : Orders.check ⟨.add, .mul, .sub, .div⟩ 9
    [⟨2, 2, 4⟩, ⟨2, 4, 8⟩, ⟨8, 2, 6⟩, ⟨6, 2, 3⟩, ⟨2, 3, 5⟩, ⟨3, 3, 9⟩, ⟨9, 2, 7⟩] = false := by
  decide

