echo "Hello!"

sleep 10

echo "I just fell asleep!"

sleep 20

ssh balloflight@s483.sureserver.com 'rm /home/balloflight/www/weather/current/*.png ; rm /home/balloflight/www/weather/current/*.txt'
