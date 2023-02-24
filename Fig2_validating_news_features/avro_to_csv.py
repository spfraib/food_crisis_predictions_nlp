#!/usr/bin/python
import csv
import sys
import argparse
import io

def genSchema(coldict):

	ss = """
{"namespace": "example.avro",
 "type": "record",
 "name": "User",
 "fields": [
	%s
 ]
}""" % ",\n".join(['{"name": "%s", "type": "string"}' % x for x in coldict])
	print(ss)
	sf = file(args.schemafile,'w')
	sf.write(ss)

csv.field_size_limit(sys.maxsize)

parser = argparse.ArgumentParser(description='Convert CSV to Avro')
parser.add_argument('filename',help='Path to CSV file');
parser.add_argument('--schemafile',help='Path to output avro schema file',default='schema.avsc');
args = parser.parse_args();

f =  file(args.filename,'r')

try:
	reader = csv.DictReader(f)
	print(reader.fieldnames)
	genSchema(reader.fieldnames)
	#for row in reader:
	#	print row['BILLADDRESS'];

finally:
	f.close()
