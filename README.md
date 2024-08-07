# Flower

# On HPC
In interactive session:
# julia +1.10.4 --project=$HOME/flower/Flower.jl --threads=1 $HOME/flower/Flower.jl/examples/main_current_folder.jl $HOME/flower/Flower.jl/examples/validation.yml

## Post-processing 
In the folder where the results are stored:

python3 -c "import plot_flower; plot_flower.plot_all_fig()" half_circle.yml

<!-- # Testing -->

