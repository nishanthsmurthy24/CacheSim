mkdir $PYTHIA_HOME/traces

# Download Traces
perl download_traces.pl --csv artifact_traces.csv --dir ../traces/

# Download PARSEC
# python3 google_drive.py 1iTb4rirR2lAPtEs-rbuT_XrlYe9cJc9L $PYTHIA_HOME/traces/
gdown 1iTb4rirR2lAPtEs-rbuT_XrlYe9cJc9L
unzip PARSEC-2.1.zip -d ../traces/
rm PARSEC-2.1.zip

# Dwonload LIGRA
# python3 google_drive.py 1FKK32PuoJYEBhEQ0shX7Ej0BgUJ3DMBq ../traces/
gdown 1FKK32PuoJYEBhEQ0shX7Ej0BgUJ3DMBq
unzip Ligra.zip -d ../traces/
rm Ligra.zip

rm -rf ../traces/__MACOSX