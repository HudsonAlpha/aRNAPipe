aRNApipe is a project-oriented pipeline for processing of RNA-seq data in high performance cluster environments. The provided framework is highly modular and has been designed to be deployen on HPC environments using IBM Platform LSF, although it can be easily migrated to any other workload manager. The main features of aRNApipe are:

- Automatization and synchronization of a broad range of RNA-seq primary analyses including quality control metrics, transcript alignment, count generation, fusion identification, and variant calling.
- Project-oriented and dynamic approach allowing users to easily update analyses to include or exclude samples or enable additional processing modules
- Use of centralized parametrization files that guarantees that all the libraries will be processed using the same workflow and the same parametrization
- Use the power of HPC clusters to distribute the workload across different nodes
- Independent setting of the computational resources assigned to each processing module to optimize the use of the available resources.
- Generation of interactive web reports for sample and project tracking
- Management of genome assemblies available to perform an analysis
