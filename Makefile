install: 
	pip install -e .; qiime dev refresh-cache

test: 
	nox -r
