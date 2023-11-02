from multiprocessing import Process,cpu_count,Pool 
from tqdm import tqdm 

class MyProcess(Process):
	def __init__(self,func=None,args_list=[],num_processes=cpu_count()):
		self.func = func
		self.args_list = args_list
		self.num_processes = num_processes
		#自己写__init__(self)会将父类的__init__覆盖，为了不丢失父类的一些属性，需要用super()加载
		super().__init__()
		

	def run_multiProcess(self):
		p = Pool(self.num_processes)
		print("总进程数: {}个".format(self.num_processes))
		result = p.map(self.func,self.args_list)
		result = list(tqdm(p.imap(self.func, self.args_list), total=len(self.args_list)))
		return result
	

def test_func(i):
	import time
	time.sleep(1)
	return  i


if __name__ == "__main__":
	import numpy as np
	import time
	
	args_list = np.arange(0,100,1)

    # 多进程
	mp_start_time = time.time()
	mp = MyProcess(func=test_func,args_list=args_list)
	result = mp.run_multiProcess()
	# print(result)
	mp_end_time = time.time()
	print("多进程耗时: {:.2f}秒".format(mp_end_time - mp_start_time))
	
	
    # 单进程
	sp_start_time = time.time()
	result = np.zeros_like(args_list)
	for i in args_list:
		result[i] = test_func(i)
		
	sp_end_time = time.time()
	print("单进程耗时: {:.2f}秒".format(sp_end_time - sp_start_time))  
	



    
    
	


