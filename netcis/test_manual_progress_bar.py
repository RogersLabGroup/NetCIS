"""
Code from ChatGPT that answered the following question:

how do I manually update the progress bar in tqdm package 
in python when I am using it with Python's multiprocessing 
imap_unordered function?
"""


import time
from multiprocessing import Pool, Manager
from tqdm import tqdm
import random

# Function to simulate work in each process
def worker(task, progress_queue):
    time.sleep(random.random())  # Simulate some work
    progress_queue.put(1)  # Put 1 in the queue to indicate completion of a task

def main():
    # List of tasks to process (replace with your actual tasks)
    tasks = range(100)
    
    # Create a Manager and Queue to share data between processes
    with Manager() as manager:
        # Queue to track progress
        progress_queue = manager.Queue()  

        # Create a tqdm progress bar
        with tqdm(total=len(tasks)) as pbar:
            
            # Create a Pool of workers
            with Pool(10) as pool:
                
                # Submit tasks asynchronously
                for task in tasks:
                    pool.apply_async(worker, args=(task, progress_queue))

                # Monitor progress and update the progress bar
                completed_tasks = 0
                while completed_tasks < len(tasks):
                    progress_queue.get()  # Wait for one task to complete
                    pbar.update(1)  # Update the progress bar
                    completed_tasks += 1
                
                pool.close()  # Close the pool to prevent new tasks
                pool.join()  # Wait for all worker processes to finish

if __name__ == "__main__":
    main()
