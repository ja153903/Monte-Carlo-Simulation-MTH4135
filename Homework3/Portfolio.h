
/*
Due to the absence of classes in C,
I opted to use structs to handle Portfolio activity
This also helped with the readability of my code
*/

struct Portfolio{
	/*
	We initialize the struct with three variables
	funds: the amount of money we have
	shares: the number of shares we hold
	mma: indicator variable denoting whether we are invested
	fully in money market account or not.
	*/
	double funds;
	double shares;
	int mma;
};

void interest(Portfolio& p, double t, double r){
	/*
	This function applies interest to our cash
	at a specific time with some rate.
	*/
	p.funds = p.funds * exp(r * t);
}

void updateShares(Portfolio& p, double price, double amount){
	/*
	This function updates the portfolio's shares and funds 
	accordingly if any changes to the investment occurs
	*/
	p.shares = p.shares + amount;
	p.funds = p.funds - price * amount;
}

void closePositions(Portfolio& p, double price){
	/*
	This function closes the portfolio's positions
	liquidating any assets held.
	*/
	p.funds = p.funds + p.shares * price;
	p.shares = 0;
}

void clearAll(Portfolio& p, double P0){
	/*
	This function resets the portfolio to its initial state
	for simulation purposes.
	*/
	p.funds = P0;
	p.shares = 0;
	p.mma = 1;
}

void strategy(Portfolio& p, double price, double t, double T, double r, char question){

	if (question == 'a'){
		interest(p, t, r);
	}

	if (question == 'b'){
		updateShares(p, price, p.funds/price);
	}

	if (question == 'c'){

		if (p.mma == 1 && price < 95){
			updateShares(p, price, p.funds/price);
			p.mma = 0;
		}

		if (p.mma == 0 && price > 105){
			updateShares(p, price, -p.shares);
			p.mma = 1;
		}
		interest(p, t, r);
	}

	if (question == 'd'){

		if (p.mma == 1 && price > 105){
			updateShares(p, price, p.funds/price);
			p.mma = 0;
		}

		if (p.mma == 0 && price < 95){
			updateShares(p, price, -p.shares);
			p.mma = 1;
		}
		interest(p, t, r);
	}

	if (question == 'e'){

		double vt = 110 * exp (-r *T);

		if (p.shares == 0 && price >= vt){
			updateShares(p, price, 1);
		}
		
		if (p.shares == 1 && price <= vt){
			updateShares(p, price, -1);
		}
		interest(p, t, r);
	}
}