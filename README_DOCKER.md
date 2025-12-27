# Docker Setup for Air Spring Calculator

This project includes Docker support for easy deployment of the API service.

## Prerequisites

- Docker installed on your system

## Building the Docker Image

```bash
docker build -t air-spring-calculator .
```

## Running the Container

```bash
docker run -d -p 5000:5000 --name air-spring-calc air-spring-calculator
```

This will:
- Start the container in detached mode
- Expose the API on port 5000
- Name the container `air-spring-calc`

## Accessing the API

Once the container is running, the API will be available at:
- `http://localhost:5000`

### Test the API:

```bash
# Health check
curl http://localhost:5000/health

# Calculate spring rate
curl -X POST http://localhost:5000/api/calculate \
  -H "Content-Type: application/json" \
  -d '{
    "travel": 160,
    "startingPressure": 56,
    "irtPressure": 86,
    "travelPoint": 80
  }'
```

## Stopping the Container

```bash
docker stop air-spring-calc
docker rm air-spring-calc
```

Or to stop and remove in one command:
```bash
docker rm -f air-spring-calc
```

## Viewing Logs

```bash
docker logs air-spring-calc
```

To follow logs in real-time:
```bash
docker logs -f air-spring-calc
```

## Environment Variables

The container uses the following environment variables:
- `FLASK_ENV`: Set to `production` (disables debug mode)
- `FLASK_APP`: Set to `api_service.py`

You can override these when running the container:
```bash
docker run -d -p 5000:5000 \
  -e FLASK_ENV=production \
  --name air-spring-calc \
  air-spring-calculator
```

## Port Configuration

The API service runs on port 5000 inside the container. To change the external port, modify the port mapping:

```bash
docker run -d -p 8080:5000 --name air-spring-calc air-spring-calculator
```

This maps external port 8080 to internal port 5000.

## Restart Policy

To automatically restart the container if it stops:
```bash
docker run -d -p 5000:5000 --restart unless-stopped --name air-spring-calc air-spring-calculator
```
