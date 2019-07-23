build:
	docker build -t napviewer .

bash:
	docker run -it -p 5001:5001 --rm --name napviewer napviewer bash

interactive:
	docker run -it -p 5001:5001 --rm -v /tmp:/NAPviewer/api/static/downloads --name napviewer napviewer /NAPviewer/run_server.sh

server:
	docker run -it -p 5001:5001 --rm -v /tmp:/NAPviewer/api/static/downloads --name napviewer napviewer /NAPviewer/run_server.sh
