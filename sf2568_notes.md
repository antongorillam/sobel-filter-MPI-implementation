# Sobel-filter MPI implementation and 

Book on sobel filter is on this  [link](https://dokumen.pub/parallel-programming-techniques-and-applications-using-networked-workstations-and-parallel-computers-2nd-ed-0131405632-9780131405639.html)

- Create project code skeleton (**Anton**).
- Parse image to matrix.
    - Either C/C++ (preferably) or python/matlab (**Anton**)
- **Implementation** Write pseudo code
    - Understand how to partition the data to each processor.
    - Understand what each processor should do.
    - How each processor should "talk" to eachother.
- **Evaluation** Come up with different experiment ideas.
    - \# processors vs time.
    - Different filters.
    - Calculate $t_{start}$ vs $t_{compute}$ + $t_{comm}$. Perhaps it's even possible to compute $t_{compute}$ and $t_{comm}$ separately.
    - Maaaybe larger filters (probably not).
    