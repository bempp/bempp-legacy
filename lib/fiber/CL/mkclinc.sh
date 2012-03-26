# Strip comments and dereference includes
# Note that we need to create an intermediate file so that xxd can
# generate its full output
#cpp -P $1 /tmp/$1

# Create the include file containing the character array and array length
# Replace 'unsigned' with 'const'
#xxd -i /tmp/$1 - | sed -e s/unsigned/const/ -e s/_tmp_//g > $1.str

xxd -i $1 - | sed -e s/unsigned/const/ > $1.str
