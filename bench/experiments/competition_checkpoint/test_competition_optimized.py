import unittest
from fragment_competition_optimized import takeover_plan
class TakeoverTests(unittest.TestCase):
 def test_target_collision_releases_old_owner(self):
  old={0:0,1:1,2:2,3:3};owners={0:0,1:1,2:1,3:1}
  p=takeover_plan(old,owners,0,{0:0,1:2})
  self.assertEqual(p['anchors'],{0:0,1:2,3:3});self.assertEqual(p['holes'],[2])
  self.assertEqual(p['eaten'],[1]);self.assertEqual(p['displaced'],[2]);self.assertEqual(old,{0:0,1:1,2:2,3:3})
 def test_prefix_protection_and_explicit_relocation(self):
  with self.assertRaises(ValueError):takeover_plan({0:0,1:1,2:2},{0:0,1:1,2:1},0,{0:1,2:2})
  p=takeover_plan({0:0,1:1,2:2},{0:0,1:1,2:1},0,{0:1,2:2},allow_b_relocation=True)
  self.assertEqual(p['anchors'],{0:1,2:2});self.assertEqual(p['holes'],[1])
 def test_noninjective_proposal_rejected(self):
  with self.assertRaises(ValueError):takeover_plan({0:0,1:1},{0:0,1:1},0,{0:0,1:0})
 def test_partially_eaten_fragment_keeps_unaffected_pairs(self):
  p=takeover_plan({i:i for i in range(6)},{0:0,1:0,2:1,3:1,4:1,5:1},0,{0:0,2:3})
  self.assertEqual(p['anchors'],{0:0,1:1,2:3,4:4,5:5});self.assertEqual(p['retained'],{4:4,5:5})
if __name__=='__main__':unittest.main()
