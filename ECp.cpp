#include "ECp.h"

CECp::CECp(void)
{
}

CECp::CECp( char * A, char* B )
{
	CoefA.clone(CZZn(A));
	CoefB.clone(CZZn(B));
}
CECp::~CECp(void)
{
}

bool CECp::operator==( CECp& b )
{
	return ((CoefA == b.CoefA)&&(CoefB == b.CoefB)) ;
}
CECPoint::CECPoint(void) : Three(3), Two(2)
{
}


CECPoint::~CECPoint(void)
{
}

CECPoint::CECPoint(CECp* c): Three(3), Two(2)
{
	curve = c;
	infinity = true;
	x.clear();
	y.clear();
}
CECPoint::CECPoint(CECp* c, char* coord_x, char* coord_y): Three(3), Two(2)
{ 
	curve = c;
	x.clone(CZZn(coord_x));
	y.clone(CZZn(coord_y));
	infinity = (x.isZero() && y.isZero());
}
CECPoint::CECPoint(CECp* c, CZZn coord_x, CZZn coord_y) : Three(3), Two(2){
	curve = c;
	x.clone(coord_x);
	y.clone(coord_y);
	infinity = (x.isZero() && y.isZero());
}
CECPoint::CECPoint(const CECPoint& p): Three(3), Two(2){
	curve = p.curve;
	infinity = p.infinity;
	x.clone(const_cast<CECPoint&>(p).x);
	y.clone(const_cast<CECPoint&>(p).y);
}
CECPoint CECPoint::operator+(const CECPoint& p2)
{
	CECPoint& p = const_cast<CECPoint&>(p2);
	if (!(*curve == *p.curve))
		return CECPoint(curve);
	if (p.infinity)
		return CECPoint(*this);
	if (infinity)
		return CECPoint(p);
	if (x == p.x) {//Either doubling or infinity
		if (!(y == p.y)) //infinity
			return CECPoint(curve);
		//Doubling
		CZZn grad( (CZZn(3) * x * x + curve->CoefA)/(CZZn(2)*y));
		CZZn rx(grad*grad - CZZn(2)*x);
		CZZn ry(y+grad*(rx-x));
		ry.negate();
		return CECPoint(curve,rx, ry);
	}
	//Adding
	CZZn grad((y-p.y)/(x-p.x));
	CZZn rx(grad*grad - x - p.x);
	CZZn ry(y + grad*(rx-x));
	ry.negate();
	return CECPoint(curve, rx, ry);
}
void CECPoint::operator+=(CECPoint& p)
{
	if ((curve != (p.curve)) && !(*curve == *p.curve))
	{
		return;
	}
	if (p.infinity)
	{
		return;
	}
	if (infinity)
	{
		x.clone(p.x);
		y.clone(p.y);
		infinity = p.infinity;
		return;
	}
	if (x == p.x) {//Either doubling or infinity
		if (!(y == p.y)) //infinity
		{
			infinity = true;
			return;
		}
		//Doubling
		CZZn grad( (Three * x * x + curve->CoefA)/(Two*y));
		CZZn rx(grad*grad - Two*x);
		CZZn ry(y+grad*(rx-x));
		ry.negate();
		x.clone(rx);
		y.clone(ry);
		infinity = false;
		return;
	}
	//Adding
	/*
	CZZn grad((y-p.y)/(x-p.x));
	CZZn rx(grad*grad - x - p.x);
	CZZn ry(y + grad*(rx-x));
	ry.negate();
	*/
	tmp_grad.clone(y);
	tmp_grad -= p.y;
	tmp_ry.clone(x);
	tmp_ry -= p.x;
	tmp_grad /= tmp_ry;
	tmp_rx.clone(tmp_grad);
	tmp_rx *= tmp_rx;
	tmp_rx -= p.x;
	tmp_rx -= x;
	tmp_ry.clone(tmp_rx);
	tmp_ry -= x;
	tmp_ry *= tmp_grad;
	x.clone(tmp_rx);
	y+= tmp_ry;
	y.negate();
	infinity = false;
	return;
}

CECPoint CECPoint::operator*(CZZ factor)
{
	CECPoint r(curve);
	for (int i=factor.getBitCount()-1;i>=0;i--) {
		r += r;
		if (factor.getBit(i) != 0)
			r += (*this);
	}
	return CECPoint(r);
}
CECPoint CECPoint::operator*(CZZn factor)
{
	CECPoint r(curve);
	for (int i=factor.getBitCount()-1;i>=0;i--) {
		r += r;
		if (factor.getBit(i) != 0)
			r += (*this);
	}
	return CECPoint(r);
}

bool CECPoint::isValid()
{
	CZZn t(y*y - x*x*x -x*curve->CoefA-curve->CoefB);
	return t.isZero();
}

bool CECPoint::operator==( const CECPoint& p )
{
	return ((infinity == true) && (p.infinity == true)) 
		|| ((infinity == false) && (p.infinity == false) && (x == p.x) && (y == p.y));
}
void CECMultiPoints::set(const CECPoint& pt, int index )
{
	points[index] = pt;
}

CECPoint& CECMultiPoints::get( int index )
{
	return points[index];
}

void CECMultiPoints::MultipleAdd( CECPoint* pts[NUM_MULTIPOINTS] )
{
	int i;
	bool diverge = false;
	for(i=0;i<NUM_MULTIPOINTS;i++)
	{
		if (pts[i]->infinity)
		{
			//Adding a null point
			diverge = true;
		}
		else if (points[i].infinity)
		{
			//I am a null point
			points[i] = *pts[i];
			diverge = true;
		}
		else if (points[i].x == pts[i]->x) {//Either doubling or infinity
			if (!(points[i].y == pts[i]->y)) //infinity
			{
				points[i].infinity = true;
				diverge = true;
			} 
			else
			{
				//Doubling
				points[i] += *pts[i];
				diverge = true;
			}
		}
		else //Adding two different points
		{
			if (diverge) // Different operations already performed on other points, so use standard add methods.
			{
				points[i] += *pts[i];
			}
		}

	}
	if (diverge) //All work has been done if diverge is true.
		return;
	
	CZZn* inverses[NUM_MULTIPOINTS];
	// Now all points additions are normal adding, use the aggregating modular inverse trick.
	for(i=0;i<NUM_MULTIPOINTS;i++)
	{
		grad[i].clone(points[i].y);
		grad[i] -= pts[i]->y;
		ry[i].clone(points[i].x);
		ry[i] -= pts[i]->x;
		if (NUM_MULTIPOINTS > 1)
			inverses[i] = &ry[i];
	}
	if (NUM_MULTIPOINTS > 1)
		CZZn::MultipleInvert(inverses, NUM_MULTIPOINTS);
	else
		grad[0] /= ry[0];
	for(i=0;i<NUM_MULTIPOINTS;i++)
	{
		if (NUM_MULTIPOINTS > 1)
			grad[i] *= *inverses[i]; //Gradient done
		rx[i].clone(grad[i]);
		rx[i] *= rx[i];
		rx[i] -= pts[i]->x;
		rx[i] -= points[i].x;
		ry[i].clone(rx[i]);
		ry[i] -= points[i].x;
		ry[i] *= grad[i];
		points[i].x.clone(rx[i]);
		points[i].y+= ry[i];
		points[i].y.negate();
		points[i].infinity = false;
	}
}

CECMultiPoints::CECMultiPoints(CECp* curve): Three(3), Two(2)
{
	this->curve = curve;
}
