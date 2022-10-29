{% extends "slurm.sh" %}

{% block header %}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
    {{- super () -}}

{% if gpus %}
#SBATCH -q gpu
#SBATCH --gres gpu:{{ gpus }}
{%- else %}
#SBATCH -q requeue
{%- endif %}

#SBATCH -N 1
#SBATCH --mail-type=ALL




echo  "Running on host" hostname
echo  "Time is" date
conda activate mosdef-study38-new

module load python/3.8

{% if gpus %}
module load cuda/10.0
nvidia-smi
{%- endif %}

{% endblock header %}

{% block body %}
    {{- super () -}}


{% endblock body %}
