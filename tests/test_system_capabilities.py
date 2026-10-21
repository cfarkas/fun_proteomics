from pathlib import Path
import sys,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from system_capabilities import recommend,GIB
class ResourceTests(unittest.TestCase):
 def test_host_and_reserve(self):
  c=recommend(128,128,1000*GIB,980*GIB)
  self.assertEqual(c['recommended_threads'],128);self.assertEqual(c['recommended_memory_gib'],882)
 def test_affinity_and_quota(self):
  c=recommend(128,16,1000*GIB,980*GIB,quotas=[4.5,8])
  self.assertEqual(c['available_cpus'],4)
 def test_cgroup_memory_headroom(self):
  c=recommend(8,8,64*GIB,60*GIB,memory_limits=[(16*GIB,6*GIB),(32*GIB,8*GIB)])
  self.assertEqual(c['max_memory_gib'],10);self.assertEqual(c['recommended_memory_gib'],9)
 def test_insufficient_memory(self):
  self.assertFalse(recommend(4,4,4*GIB,3*GIB)['supported'])
