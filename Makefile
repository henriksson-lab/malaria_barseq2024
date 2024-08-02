docker:
	docker build -t malaria_barseq2024 .
docker_testrun:
	docker run --rm -p 3838:3838 malaria_barseq2024
docker_push:
        #docker images
	#docker tag 045a4cba1107 mahogny83/malaria_barseq2024:20240416-095800
	#docker image push mahogny83/malaria_barseq2024:20240416-095800
	#need to update tag
