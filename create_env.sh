conda env create -f conda-env.yml python=3.7 -p ~/env/$ENV_NAME

# dump env
# conda env export --from-history | sed "/^prefix:/d" > conda-env.yml
